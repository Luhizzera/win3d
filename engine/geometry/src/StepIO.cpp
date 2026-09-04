#include "aether/geometry/StepIO.hpp"

#include "aether/mesh/PolygonTriangulation2D.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

#ifdef AETHER_HAVE_OPENCASCADE
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IGESControl_Reader.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Reader.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#endif

namespace aether::geometry {

bool stepIoHasOpenCascade() {
#ifdef AETHER_HAVE_OPENCASCADE
    return true;
#else
    return false;
#endif
}

#ifdef AETHER_HAVE_OPENCASCADE

using core::Vector3;

namespace {

// Tessellates every face of `shape` (BRepMesh_IncrementalMesh handles a
// planar and a curved face identically -- the whole reason this path
// exists over the tokenizer below) and packs the result into a
// TriangleMesh. A face's own TopLoc_Location carries whatever assembly-
// level placement was applied above it; applying that transform to each
// node is what makes the triangulation correct in the shape's own
// coordinate system rather than the face's local one. TopAbs_REVERSED
// flips a face's sense relative to its underlying surface's natural
// normal, so triangle winding is flipped here to keep every triangle's
// normal pointing outward -- the curved-geometry equivalent of the
// tokenizer's own care for FACE_OUTER_BOUND's .T./.F. sense flag.
StepLoadResult tessellateShape(const TopoDS_Shape& shape) {
    StepLoadResult result;
    if (shape.IsNull()) {
        result.unsupportedFeatures.push_back("file transferred no shape");
        return result;
    }

    BRepMesh_IncrementalMesh mesh(shape, /*theLinDeflection=*/1e-3, /*isRelative=*/true,
                                  /*theAngDeflection=*/0.3);
    (void)mesh; // constructor performs the meshing; the object itself is not queried again

    for (TopExp_Explorer explorer(shape, TopAbs_FACE); explorer.More(); explorer.Next()) {
        const TopoDS_Face& face = TopoDS::Face(explorer.Current());
        TopLoc_Location location;
        const occ::handle<Poly_Triangulation>& triangulation = BRep_Tool::Triangulation(face, location);
        if (triangulation.IsNull()) {
            result.unsupportedFeatures.push_back("a face produced no triangulation");
            continue;
        }
        const gp_Trsf& transform = location.Transformation();
        const bool reversed = face.Orientation() == TopAbs_REVERSED;

        const std::size_t firstVertex = result.mesh.vertexCount();
        for (int i = 1; i <= triangulation->NbNodes(); ++i) {
            const gp_Pnt p = triangulation->Node(i).Transformed(transform);
            result.mesh.addVertex(Vector3(p.X(), p.Y(), p.Z()));
        }
        for (int i = 1; i <= triangulation->NbTriangles(); ++i) {
            int n1 = 0, n2 = 0, n3 = 0;
            triangulation->Triangle(i).Get(n1, n2, n3);
            const std::size_t a = firstVertex + static_cast<std::size_t>(n1 - 1);
            const std::size_t b = firstVertex + static_cast<std::size_t>(n2 - 1);
            const std::size_t c = firstVertex + static_cast<std::size_t>(n3 - 1);
            if (reversed) {
                result.mesh.addTriangle(a, c, b);
            } else {
                result.mesh.addTriangle(a, b, c);
            }
        }
    }

    if (result.mesh.triangleCount() == 0) {
        result.unsupportedFeatures.push_back("no triangulated face found in this file's shape");
    }
    return result;
}

} // namespace

StepLoadResult loadStep(const std::string& path) {
    STEPControl_Reader reader;
    if (reader.ReadFile(path.c_str()) != IFSelect_RetDone) {
        throw std::runtime_error("loadStep: cannot open or parse '" + path + "' as a STEP file");
    }
    reader.TransferRoots();
    return tessellateShape(reader.OneShape());
}

StepLoadResult loadIges(const std::string& path) {
    IGESControl_Reader reader;
    if (reader.ReadFile(path.c_str()) != IFSelect_RetDone) {
        throw std::runtime_error("loadIges: cannot open or parse '" + path + "' as an IGES file");
    }
    reader.TransferRoots();
    return tessellateShape(reader.OneShape());
}

#else // !AETHER_HAVE_OPENCASCADE

namespace {

using core::Vector3;

// -- Part 21 (ISO 10303-21) tokenizer -------------------------------------
//
// Purely syntactic: this section knows the *text format*'s grammar (quoted
// strings with doubled-quote escaping, parenthesised lists, entity
// references, simple vs. complex entity instances) and nothing about any
// EXPRESS schema. That split matters -- the text format is small, stable
// (unchanged since 1994) and fully specified, so it can be implemented and
// tested with confidence; the schema is enormous, and only the specific
// entities this loader actually interprets (below) were individually
// checked against the published standard before being relied on.

struct StepValue {
    enum class Kind { String, Reference, List, Other };
    Kind kind = Kind::Other;
    std::string text;             // String: unescaped text. Reference: digits only. Other: raw token.
    std::vector<StepValue> items; // List only.
};

struct StepEntityInstance {
    // More than one entry means a Part 21 *complex* entity instance
    // (`#1=(TYPE_A(...)TYPE_B(...));`), the syntax AP242's tessellated
    // representations use and this loader deliberately does not interpret
    // -- see this file's header comment. Its type names are still recorded
    // so the loader can name it in a diagnostic rather than silently
    // ignore it.
    std::vector<std::string> typeNames;
    std::vector<std::vector<StepValue>> argsByType;
};

class Part21Parser {
public:
    explicit Part21Parser(const std::string& text) : text_(text) {}

    std::map<long long, StepEntityInstance> parse() {
        std::map<long long, StepEntityInstance> entities;
        while (true) {
            skipToNextEntityMarker();
            if (atEnd()) {
                break;
            }
            const long long id = parseEntityId();
            skipWhitespaceAndComments();
            expect('=');
            skipWhitespaceAndComments();
            StepEntityInstance entity = parseEntityInstance();
            skipWhitespaceAndComments();
            expect(';');
            entities[id] = std::move(entity);
        }
        return entities;
    }

private:
    // Advances past anything that is not the start of an `#<digits>`
    // reference -- comments, header macros, whitespace -- while still
    // respecting quoted strings, so a literal '#' inside a filename or
    // description string is never mistaken for an entity marker. Leaves
    // the cursor exactly on the '#' when one is found.
    void skipToNextEntityMarker() {
        while (!atEnd()) {
            skipWhitespaceAndComments();
            if (atEnd()) {
                return;
            }
            if (peek() == '\'') {
                skipQuotedString();
                continue;
            }
            if (peek() == '#' && pos_ + 1 < text_.size() && std::isdigit(static_cast<unsigned char>(text_[pos_ + 1]))) {
                return;
            }
            ++pos_;
        }
    }

    long long parseEntityId() {
        expect('#');
        const std::size_t start = pos_;
        while (!atEnd() && std::isdigit(static_cast<unsigned char>(peek()))) {
            ++pos_;
        }
        return std::stoll(text_.substr(start, pos_ - start));
    }

    // Identifier := TYPE(args) | ( Identifier(args) Identifier(args) ... )
    StepEntityInstance parseEntityInstance() {
        StepEntityInstance entity;
        if (peek() == '(') {
            ++pos_; // consume the complex-instance's own outer '('
            skipWhitespaceAndComments();
            while (peek() != ')') {
                entity.typeNames.push_back(parseIdentifier());
                skipWhitespaceAndComments();
                entity.argsByType.push_back(parseArgList());
                skipWhitespaceAndComments();
            }
            ++pos_; // consume the outer ')'
        } else {
            entity.typeNames.push_back(parseIdentifier());
            skipWhitespaceAndComments();
            entity.argsByType.push_back(parseArgList());
        }
        return entity;
    }

    std::string parseIdentifier() {
        const std::size_t start = pos_;
        while (!atEnd() && (std::isalnum(static_cast<unsigned char>(peek())) || peek() == '_')) {
            ++pos_;
        }
        if (pos_ == start) {
            throw std::runtime_error("loadStep: expected an entity type name at offset " +
                                      std::to_string(pos_));
        }
        return text_.substr(start, pos_ - start);
    }

    // Consumes '(' value (',' value)* ')'. An empty list '()' is valid.
    std::vector<StepValue> parseArgList() {
        expect('(');
        std::vector<StepValue> values;
        skipWhitespaceAndComments();
        if (peek() == ')') {
            ++pos_;
            return values;
        }
        while (true) {
            values.push_back(parseValue());
            skipWhitespaceAndComments();
            if (peek() == ',') {
                ++pos_;
                skipWhitespaceAndComments();
                continue;
            }
            expect(')');
            break;
        }
        return values;
    }

    StepValue parseValue() {
        skipWhitespaceAndComments();
        StepValue value;
        if (peek() == '\'') {
            value.kind = StepValue::Kind::String;
            value.text = readQuotedString();
        } else if (peek() == '(') {
            value.kind = StepValue::Kind::List;
            value.items = parseArgList();
        } else if (peek() == '#') {
            value.kind = StepValue::Kind::Reference;
            ++pos_;
            const std::size_t start = pos_;
            while (!atEnd() && std::isdigit(static_cast<unsigned char>(peek()))) {
                ++pos_;
            }
            value.text = text_.substr(start, pos_ - start);
        } else {
            // A bare token: number, enumeration (.T./.F./.PLANE.-style),
            // '$' (unset) or '*' (derived) -- all of them run until the
            // next ',' or ')' at this level, since none can contain one.
            const std::size_t start = pos_;
            while (!atEnd() && peek() != ',' && peek() != ')') {
                ++pos_;
            }
            value.kind = StepValue::Kind::Other;
            value.text = text_.substr(start, pos_ - start);
            // Trim trailing whitespace the scan above may have absorbed
            // before a ',' or ')' that had space in front of it.
            while (!value.text.empty() && std::isspace(static_cast<unsigned char>(value.text.back()))) {
                value.text.pop_back();
            }
        }
        return value;
    }

    void skipQuotedString() { readQuotedString(); }

    // STEP escapes an embedded quote by doubling it (`'it''s'` -> `it's`).
    std::string readQuotedString() {
        expect('\'');
        std::string result;
        while (true) {
            if (atEnd()) {
                throw std::runtime_error("loadStep: unterminated string literal");
            }
            const char c = text_[pos_++];
            if (c != '\'') {
                result.push_back(c);
                continue;
            }
            if (!atEnd() && peek() == '\'') {
                result.push_back('\'');
                ++pos_;
                continue;
            }
            break; // the closing quote
        }
        return result;
    }

    void skipWhitespaceAndComments() {
        while (!atEnd()) {
            if (std::isspace(static_cast<unsigned char>(peek()))) {
                ++pos_;
            } else if (peek() == '/' && pos_ + 1 < text_.size() && text_[pos_ + 1] == '*') {
                pos_ += 2;
                while (pos_ + 1 < text_.size() && !(text_[pos_] == '*' && text_[pos_ + 1] == '/')) {
                    ++pos_;
                }
                pos_ = std::min(pos_ + 2, text_.size());
            } else {
                break;
            }
        }
    }

    bool atEnd() const { return pos_ >= text_.size(); }
    char peek() const { return text_[pos_]; }
    void expect(char c) {
        if (atEnd() || peek() != c) {
            throw std::runtime_error(std::string("loadStep: expected '") + c + "' at offset " +
                                      std::to_string(pos_));
        }
        ++pos_;
    }

    const std::string& text_;
    std::size_t pos_ = 0;
};

// -- Schema-level extraction: the faceted-BREP path only ------------------
//
// Every entity shape below (attribute names, types, order) was checked
// against ISO 10303-42's published EXPRESS text before being relied on --
// see StepIO.hpp's header comment for which ones and why this stops here.

const StepValue* asList(const StepValue& v) {
    return v.kind == StepValue::Kind::List ? &v : nullptr;
}
bool asReference(const StepValue& v, long long& out) {
    if (v.kind != StepValue::Kind::Reference) {
        return false;
    }
    out = std::stoll(v.text);
    return true;
}
bool asBool(const StepValue& v, bool& out) {
    if (v.kind != StepValue::Kind::Other) {
        return false;
    }
    if (v.text == ".T.") {
        out = true;
        return true;
    }
    if (v.text == ".F.") {
        out = false;
        return true;
    }
    return false;
}
bool asDouble(const StepValue& v, double& out) {
    if (v.kind != StepValue::Kind::Other || v.text.empty()) {
        return false;
    }
    try {
        std::size_t consumed = 0;
        out = std::stod(v.text, &consumed);
        return consumed == v.text.size();
    } catch (...) {
        return false;
    }
}

// Type names that mean "this file contains geometry only a real CAD kernel
// can evaluate", collected once so the surface scan and the per-face
// lookup failures below can both point to the same explanation. Not
// exhaustive by design -- it names the common cases, and any entity that
// simply fails to resolve as a POLY_LOOP is reported on its own terms
// regardless of whether its actual type is on this list.
bool namesCurvedOrTessellatedGeometry(const std::string& typeName) {
    static const std::vector<std::string> kMarkers = {
        "B_SPLINE_SURFACE", "B_SPLINE_SURFACE_WITH_KNOTS", "RATIONAL_B_SPLINE_SURFACE",
        "BOUNDED_SURFACE", "SURFACE_OF_REVOLUTION", "SURFACE_OF_LINEAR_EXTRUSION",
        "CYLINDRICAL_SURFACE", "CONICAL_SURFACE", "SPHERICAL_SURFACE", "TOROIDAL_SURFACE",
        "ADVANCED_FACE", "EDGE_LOOP",
        "TESSELLATED_SURFACE_SET", "TRIANGULATED_SURFACE_SET", "COMPLEX_TRIANGULATED_SURFACE_SET",
        "TESSELLATED_SHELL", "TESSELLATED_SOLID", "COORDINATES_LIST",
    };
    for (const std::string& marker : kMarkers) {
        if (typeName == marker) {
            return true;
        }
    }
    return false;
}

} // namespace

StepLoadResult loadStep(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("loadStep: cannot open '" + path + "'");
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    const std::string text = buffer.str();

    if (text.find("ISO-10303-21") == std::string::npos) {
        throw std::runtime_error("loadStep: '" + path +
                                  "' does not look like an ISO-10303-21 (STEP) file "
                                  "(missing the ISO-10303-21 header marker)");
    }

    Part21Parser parser(text);
    const std::map<long long, StepEntityInstance> entities = parser.parse();

    StepLoadResult result;

    // A quick, independent pass over every entity's type name(s): even a
    // face this loader never reaches (because the walk below stopped
    // earlier) is worth naming, so the diagnostic reflects the whole file
    // rather than only the part that was traversed.
    for (const auto& [id, entity] : entities) {
        for (const std::string& typeName : entity.typeNames) {
            if (namesCurvedOrTessellatedGeometry(typeName)) {
                result.unsupportedFeatures.push_back("#" + std::to_string(id) + ": " + typeName +
                                                      " (needs a full CAD kernel, not a faceted "
                                                      "representation this loader can interpret)");
            }
        }
    }

    // -- CARTESIAN_POINT('name', (x, y, z)) -> mesh vertex, lazily --------
    // Lazy rather than eager: only points a POLY_LOOP actually references
    // become mesh vertices, so an unrelated construction point elsewhere
    // in the file (there are usually many -- placements, axes, curve
    // control points) never shows up as an orphaned, disconnected vertex.
    std::map<long long, std::size_t> pointToVertex;
    const auto resolvePoint = [&](long long id) -> bool {
        if (pointToVertex.count(id) != 0) {
            return true;
        }
        const auto it = entities.find(id);
        if (it == entities.end() || it->second.typeNames.size() != 1 ||
            it->second.typeNames[0] != "CARTESIAN_POINT") {
            return false;
        }
        const auto& args = it->second.argsByType[0];
        const StepValue* coords;
        if (args.size() != 2 || (coords = asList(args[1])) == nullptr || coords->items.size() != 3) {
            return false;
        }
        double x = 0.0, y = 0.0, z = 0.0;
        if (!asDouble(coords->items[0], x) || !asDouble(coords->items[1], y) ||
            !asDouble(coords->items[2], z)) {
            return false;
        }
        pointToVertex[id] = result.mesh.addVertex(Vector3(x, y, z));
        return true;
    };

    // -- POLY_LOOP('name', (#p1, #p2, ...)) -> ordered mesh-vertex list ---
    const auto resolveLoop = [&](long long id, std::vector<std::size_t>& out) -> bool {
        const auto it = entities.find(id);
        if (it == entities.end() || it->second.typeNames.size() != 1 ||
            it->second.typeNames[0] != "POLY_LOOP") {
            return false; // most commonly: it is an EDGE_LOOP -- a curved boundary
        }
        const auto& args = it->second.argsByType[0];
        const StepValue* polygon;
        if (args.size() != 2 || (polygon = asList(args[1])) == nullptr || polygon->items.size() < 3) {
            return false;
        }
        out.clear();
        out.reserve(polygon->items.size());
        for (const StepValue& item : polygon->items) {
            long long pointId = 0;
            if (!asReference(item, pointId) || !resolvePoint(pointId)) {
                return false;
            }
            out.push_back(pointToVertex[pointId]);
        }
        return true;
    };

    // Triangulates one planar face bound and appends its triangles.
    // `loopVertices` is already oriented (FACE_BOUND's own orientation
    // flag has been applied by the caller): Newell's method gets this
    // face's normal from the *actual* (possibly slightly non-planar, in
    // real files) point set rather than a single 3-point cross product,
    // which is unreliable if those three happen to be near-collinear or
    // the loop is concave there.
    const auto triangulateFace = [&](const std::vector<std::size_t>& loopVertices) {
        Vector3 normal;
        const std::size_t n = loopVertices.size();
        for (std::size_t i = 0; i < n; ++i) {
            const Vector3& a = result.mesh.vertex(loopVertices[i]);
            const Vector3& b = result.mesh.vertex(loopVertices[(i + 1) % n]);
            normal.x += (a.y - b.y) * (a.z + b.z);
            normal.y += (a.z - b.z) * (a.x + b.x);
            normal.z += (a.x - b.x) * (a.y + b.y);
        }
        if (normal.normSquared() == 0.0) {
            return; // degenerate (collinear) loop: nothing to triangulate
        }
        normal = normal.normalized();
        // An arbitrary in-plane basis: project world-X onto the plane
        // unless the normal is too close to world-X, then use world-Y --
        // the standard "pick whichever axis isn't nearly parallel" trick.
        Vector3 reference = std::fabs(normal.x) < 0.9 ? Vector3(1.0, 0.0, 0.0) : Vector3(0.0, 1.0, 0.0);
        const Vector3 u = (reference - normal * normal.dot(reference)).normalized();
        const Vector3 v = normal.cross(u);

        mesh::PolygonTriangulation2D polygon2D;
        const Vector3& origin = result.mesh.vertex(loopVertices[0]);
        for (const std::size_t vertexIndex : loopVertices) {
            const Vector3 offset = result.mesh.vertex(vertexIndex) - origin;
            polygon2D.addVertex(offset.dot(u), offset.dot(v));
        }
        // PolygonTriangulation2D requires CCW winding; Newell's normal is
        // already the outward one implied by loopVertices' own order, so a
        // negative signed area means (u, v) came out left-handed for this
        // particular loop and needs reversing rather than the loop itself.
        if (polygon2D.polygonArea() < 0.0) {
            mesh::PolygonTriangulation2D reversed;
            for (auto it = loopVertices.rbegin(); it != loopVertices.rend(); ++it) {
                const Vector3 offset = result.mesh.vertex(*it) - origin;
                reversed.addVertex(offset.dot(u), offset.dot(v));
            }
            reversed.triangulate();
            for (std::size_t t = 0; t < reversed.triangleCount(); ++t) {
                const auto& tri = reversed.triangle(t);
                const std::size_t a = loopVertices[loopVertices.size() - 1 - tri.vertices[0]];
                const std::size_t b = loopVertices[loopVertices.size() - 1 - tri.vertices[1]];
                const std::size_t c = loopVertices[loopVertices.size() - 1 - tri.vertices[2]];
                result.mesh.addTriangle(a, b, c);
            }
            return;
        }
        polygon2D.triangulate();
        for (std::size_t t = 0; t < polygon2D.triangleCount(); ++t) {
            const auto& tri = polygon2D.triangle(t);
            result.mesh.addTriangle(loopVertices[tri.vertices[0]], loopVertices[tri.vertices[1]],
                                     loopVertices[tri.vertices[2]]);
        }
    };

    // -- FACE_BOUND / FACE_OUTER_BOUND('name', #loop, orientation) --------
    // FACE_OUTER_BOUND adds no attributes of its own (verified: it is
    // declared with an empty body, `SUBTYPE OF (face_bound); END_ENTITY;`),
    // so both share this same 3-argument shape.
    const auto resolveFaceBound = [&](long long id, std::vector<std::size_t>& outLoop) -> bool {
        const auto it = entities.find(id);
        if (it == entities.end() || it->second.typeNames.size() != 1) {
            return false;
        }
        const std::string& type = it->second.typeNames[0];
        if (type != "FACE_BOUND" && type != "FACE_OUTER_BOUND") {
            return false;
        }
        const auto& args = it->second.argsByType[0];
        long long loopId = 0;
        bool orientation = true;
        if (args.size() != 3 || !asReference(args[1], loopId) || !asBool(args[2], orientation)) {
            return false;
        }
        if (!resolveLoop(loopId, outLoop)) {
            return false;
        }
        if (!orientation) {
            std::reverse(outLoop.begin(), outLoop.end());
        }
        return true;
    };

    // -- Walk every FACETED_BREP / MANIFOLD_SOLID_BREP in the file --------
    // Collected directly by type name rather than navigated from a
    // product's own shape representation: correct for the overwhelmingly
    // common single-part export, and it means a multi-solid file loads the
    // union of every solid's own geometry without applying any assembly
    // placement transform between them -- out of scope for this loader,
    // and worth stating rather than silently assuming every file is a
    // single part in one coordinate system.
    std::size_t solidsFound = 0;
    for (const auto& [id, entity] : entities) {
        if (entity.typeNames.size() != 1) {
            continue;
        }
        const std::string& type = entity.typeNames[0];
        if (type != "MANIFOLD_SOLID_BREP" && type != "FACETED_BREP") {
            continue;
        }
        ++solidsFound;
        const auto& args = entity.argsByType[0];
        long long shellId = 0;
        if (args.size() != 2 || !asReference(args[1], shellId)) {
            result.unsupportedFeatures.push_back("#" + std::to_string(id) + ": " + type +
                                                  " (malformed 'outer' reference)");
            continue;
        }
        const auto shellIt = entities.find(shellId);
        if (shellIt == entities.end() || shellIt->second.typeNames.size() != 1 ||
            shellIt->second.typeNames[0] != "CLOSED_SHELL") {
            result.unsupportedFeatures.push_back("#" + std::to_string(id) + ": " + type +
                                                  " (its outer shell is not a plain CLOSED_SHELL)");
            continue;
        }
        const StepValue* facesList;
        if (shellIt->second.argsByType[0].size() != 2 ||
            (facesList = asList(shellIt->second.argsByType[0][1])) == nullptr) {
            result.unsupportedFeatures.push_back("#" + std::to_string(shellId) +
                                                  ": CLOSED_SHELL (malformed face set)");
            continue;
        }

        for (const StepValue& faceRef : facesList->items) {
            long long faceId = 0;
            if (!asReference(faceRef, faceId)) {
                continue;
            }
            const auto faceIt = entities.find(faceId);
            if (faceIt == entities.end() || faceIt->second.typeNames.size() != 1 ||
                faceIt->second.typeNames[0] != "FACE") {
                result.unsupportedFeatures.push_back(
                    "#" + std::to_string(faceId) + ": not a plain FACE (likely ADVANCED_FACE with "
                    "curved geometry)");
                continue;
            }
            const StepValue* boundsList;
            if (faceIt->second.argsByType[0].size() != 2 ||
                (boundsList = asList(faceIt->second.argsByType[0][1])) == nullptr ||
                boundsList->items.empty()) {
                result.unsupportedFeatures.push_back("#" + std::to_string(faceId) +
                                                      ": FACE (malformed bound set)");
                continue;
            }
            if (boundsList->items.size() > 1) {
                result.unsupportedFeatures.push_back(
                    "#" + std::to_string(faceId) +
                    ": FACE has more than one bound (a hole) -- not supported, "
                    "PolygonTriangulation2D does not triangulate polygons with holes");
                continue;
            }
            long long boundId = 0;
            std::vector<std::size_t> loopVertices;
            if (!asReference(boundsList->items[0], boundId) || !resolveFaceBound(boundId, loopVertices)) {
                result.unsupportedFeatures.push_back(
                    "#" + std::to_string(faceId) +
                    ": face bound does not resolve to a planar POLY_LOOP (needs a CAD kernel)");
                continue;
            }
            triangulateFace(loopVertices);
        }
    }

    if (solidsFound == 0 && result.mesh.triangleCount() == 0 && result.unsupportedFeatures.empty()) {
        result.unsupportedFeatures.push_back(
            "no MANIFOLD_SOLID_BREP or FACETED_BREP entity found in this file");
    }

    return result;
}

StepLoadResult loadIges(const std::string&) {
    // No fallback exists for IGES the way the tokenizer above is one for
    // STEP's faceted subset -- IGES's entity format is different enough
    // from Part 21 that nothing here would be reusable, so this build
    // configuration has no partial capability to offer. An honest runtime
    // error, not a silently empty mesh or a missing symbol at link time
    // (this function always exists -- see StepIO.hpp's own comment on why
    // AETHER_HAVE_OPENCASCADE never appears there).
    throw std::runtime_error(
        "loadIges: built without OpenCASCADE support -- configure with "
        "-DCMAKE_TOOLCHAIN_FILE=<vcpkg-clone>/scripts/buildsystems/vcpkg.cmake to enable it");
}

#endif // AETHER_HAVE_OPENCASCADE

} // namespace aether::geometry
