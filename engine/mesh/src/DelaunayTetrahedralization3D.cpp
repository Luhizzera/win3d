#include "aether/mesh/DelaunayTetrahedralization3D.hpp"

#include "aether/mesh/RobustPredicates.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>

namespace aether::mesh {

using aether::core::Vector3;

namespace {

// Makes a positively-oriented tetrahedron (signed volume > 0) from four
// point indices, given the actual point coordinates. Orientation is
// decided by orientation3D()'s exact sign, not a plain double computation
// -- see RobustPredicates.hpp for why (this exact test is what root-fixes
// the class's documented near-degenerate-tetrahedron bug, alongside the
// in-sphere test inCircumsphere() now also uses).
DelaunayTetrahedralization3D::Tetrahedron makeTetrahedron(std::size_t a, std::size_t b, std::size_t c,
                                                           std::size_t d,
                                                           const std::vector<Vector3>& points) {
    const int sign = orientation3D(points[a], points[b], points[c], points[d]);
    if (sign < 0) {
        return {{a, b, d, c}}; // swap two vertices to flip the sign
    }
    return {{a, b, c, d}};
}

bool tetrahedronHasVertex(const DelaunayTetrahedralization3D::Tetrahedron& t, std::size_t v) {
    return t.vertices[0] == v || t.vertices[1] == v || t.vertices[2] == v || t.vertices[3] == v;
}

struct Face {
    std::size_t a;
    std::size_t b;
    std::size_t c;
};

// True if the two faces share the same 3 vertices, regardless of order.
bool sameFace(const Face& lhs, const Face& rhs) {
    const std::size_t l[3] = {lhs.a, lhs.b, lhs.c};
    const std::size_t r[3] = {rhs.a, rhs.b, rhs.c};
    for (std::size_t li : l) {
        bool found = false;
        for (std::size_t ri : r) {
            if (li == ri) {
                found = true;
                break;
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}

} // namespace

std::size_t DelaunayTetrahedralization3D::addPoint(double x, double y, double z) {
    points_.emplace_back(x, y, z);
    return points_.size() - 1;
}

namespace {

// The 4 faces of a tetrahedron, indexed so faces[li] is the face opposite
// local vertex li -- same convention makeTetrahedron()'s winding relies on.
std::array<Face, 4> tetrahedronFaces(const DelaunayTetrahedralization3D::Tetrahedron& t) {
    const auto& v = t.vertices;
    return {{{v[1], v[2], v[3]}, {v[0], v[2], v[3]}, {v[0], v[1], v[3]}, {v[0], v[1], v[2]}}};
}

// Searches every tetrahedron except `exclude` for one owning `face` (any
// vertex order). Returns its index, or -1 if `face` is a genuine hull face
// (owned by only one tetrahedron in the whole mesh).
long long findOtherFaceOwner(const std::vector<DelaunayTetrahedralization3D::Tetrahedron>& tetrahedra,
                              const Face& face, std::size_t exclude) {
    for (std::size_t t = 0; t < tetrahedra.size(); ++t) {
        if (t == exclude) {
            continue;
        }
        for (const Face& of : tetrahedronFaces(tetrahedra[t])) {
            if (sameFace(face, of)) {
                return static_cast<long long>(t);
            }
        }
    }
    return -1;
}

} // namespace

namespace {

// True iff `face` (the face of `owner` opposite its local vertex `li`) is
// visible from p -- p is on the outer side of the face's plane. See the
// class header / findCavity()'s own comment for the derivation (verified
// numerically before implementing): substituting p for the opposite
// vertex in the tetrahedron's own positively-oriented vertex order flips
// the signed volume negative exactly when p has crossed to the outside.
bool isFaceVisible(const DelaunayTetrahedralization3D::Tetrahedron& owner, int li, std::size_t p,
                    const std::vector<Vector3>& points) {
    std::array<std::size_t, 4> substituted = owner.vertices;
    substituted[li] = p;
    return orientation3D(points[substituted[0]], points[substituted[1]], points[substituted[2]],
                          points[substituted[3]]) < 0;
}

// Canonical (order-independent) key for a face, used by missingFacets(),
// hasFace() and removeRegion()'s wall lookup.
std::array<std::size_t, 3> sortedFace(const std::array<std::size_t, 3>& f) {
    std::array<std::size_t, 3> s = f;
    std::sort(s.begin(), s.end());
    return s;
}

// The tetrahedron owning `face` (any winding) together with its one vertex
// not in `face` -- the "apex" over that face, needed by
// tryFlipCoplanarQuadDiagonal() to identify the two tetrahedra straddling a
// coplanar quad. tetIndex is -1 if no tetrahedron owns the face.
struct FaceOwner {
    long long tetIndex = -1;
    std::size_t apex = 0;
};

FaceOwner findFaceOwnerWithApex(const std::vector<DelaunayTetrahedralization3D::Tetrahedron>& tetrahedra,
                                 const std::array<std::size_t, 3>& face) {
    const std::array<std::size_t, 3> key = sortedFace(face);
    for (std::size_t t = 0; t < tetrahedra.size(); ++t) {
        const auto faces = tetrahedronFaces(tetrahedra[t]);
        for (int li = 0; li < 4; ++li) {
            if (sortedFace({faces[li].a, faces[li].b, faces[li].c}) == key) {
                return {static_cast<long long>(t), tetrahedra[t].vertices[li]};
            }
        }
    }
    return {};
}

// True iff p lies inside or on tetrahedron t (which is positively oriented
// by construction, same requirement as every other predicate here). Same
// vertex-substitution trick as isFaceVisible(), just checking the sign
// stays non-negative at all four faces simultaneously instead of negative
// at one.
bool pointInTetrahedron(const DelaunayTetrahedralization3D::Tetrahedron& t, const Vector3& p,
                         const std::vector<Vector3>& points) {
    const auto& v = t.vertices;
    if (orientation3D(p, points[v[1]], points[v[2]], points[v[3]]) < 0) {
        return false;
    }
    if (orientation3D(points[v[0]], p, points[v[2]], points[v[3]]) < 0) {
        return false;
    }
    if (orientation3D(points[v[0]], points[v[1]], p, points[v[3]]) < 0) {
        return false;
    }
    if (orientation3D(points[v[0]], points[v[1]], points[v[2]], p) < 0) {
        return false;
    }
    return true;
}

} // namespace

void DelaunayTetrahedralization3D::findCavity(std::size_t p, bool considerHullVisibility,
                                               std::vector<std::size_t>& badIndices,
                                               std::vector<std::array<std::size_t, 3>>& boundary) const {
    badIndices.clear();
    boundary.clear();

    // Two independent kinds of "invalidated by p", combined below:
    //
    // (1) An ordinary REAL tetrahedron is bad iff p lies inside its
    //     circumsphere -- the plain Delaunay criterion. This alone
    //     already covers every point *inside* the current hull.
    //
    // (2) Every hull face (owned by only one real tetrahedron) has an
    //     implicit "ghost" neighbor on its far side, standing in for the
    //     unbounded region beyond the hull -- the standard trick real
    //     Delaunay implementations use for an infinite/absent neighbor. A
    //     ghost is bad iff p is visible through its face (p on the outer
    //     side), the natural limit of "p inside its circumsphere" as that
    //     far vertex recedes to infinity. This is what makes points
    //     *outside* the hull work.
    //
    // Crucially, a hull face's real owner being good does NOT make its
    // ghost good too, and vice versa -- they are judged independently.
    // (An earlier version of this function conflated the two, marking the
    // whole real tetrahedron bad whenever its ghost was bad; caught by
    // cross-checking a hand-derived example against scipy's qhull, which
    // showed the real tetrahedron must stay untouched when p is outside
    // its circumsphere, with only a single new cap tetrahedron added --
    // not a full re-fan of all its non-visible faces.)
    std::vector<bool> badReal(tetrahedra_.size(), false);
    for (std::size_t t = 0; t < tetrahedra_.size(); ++t) {
        if (inCircumsphere(tetrahedra_[t], p) > 0.0) {
            badReal[t] = true;
        }
    }

    for (std::size_t t = 0; t < tetrahedra_.size(); ++t) {
        if (badReal[t]) {
            badIndices.push_back(t);
        }
    }

    // Faces of bad real tetrahedra: boundary iff the other side is good
    // (a good real neighbor, or -- when considerHullVisibility is set -- a
    // good, non-visible ghost; otherwise every hull face is unconditionally
    // boundary, tetrahedralize()'s original behavior, always correct there
    // since every point it inserts is strictly interior to the enclosing
    // super-tetrahedron).
    for (std::size_t bi : badIndices) {
        for (int li = 0; li < 4; ++li) {
            const Face f = tetrahedronFaces(tetrahedra_[bi])[li];
            const long long other = findOtherFaceOwner(tetrahedra_, f, bi);
            if (other == -1) {
                if (!considerHullVisibility || !isFaceVisible(tetrahedra_[bi], li, p, points_)) {
                    boundary.push_back({f.a, f.b, f.c});
                }
                // else: real bad, ghost also bad -> discarded.
            } else if (!badReal[static_cast<std::size_t>(other)]) {
                boundary.push_back({f.a, f.b, f.c});
            }
            // else: shared with another bad real tetrahedron, interior to
            // the cavity, discarded.
        }
    }

    if (!considerHullVisibility) {
        return;
    }

    // Hull faces whose ghost is bad (visible) but whose real owner is
    // good: the owner stays untouched (not in badIndices), and a single
    // new cap tetrahedron (face + p) is added alongside it.
    for (std::size_t t = 0; t < tetrahedra_.size(); ++t) {
        if (badReal[t]) {
            continue; // already handled above
        }
        for (int li = 0; li < 4; ++li) {
            const Face f = tetrahedronFaces(tetrahedra_[t])[li];
            if (findOtherFaceOwner(tetrahedra_, f, t) != -1) {
                continue; // not a hull face
            }
            if (isFaceVisible(tetrahedra_[t], li, p, points_)) {
                boundary.push_back({f.a, f.b, f.c});
            }
        }
    }
}

// Routed through inSphere3D's exact arithmetic, not a plain double 4x4
// determinant -- this is the specific predicate root-caused (see the
// class header comment) as the source of a real, if rare, gap bug: a
// near-degenerate ("sliver") input tetrahedron made the old plain-double
// determinant computation unreliable. Still returns a double (exactly
// -1.0, 0.0, or 1.0) rather than changing the return type, so every
// existing `> 0.0` / `> tolerance` call site keeps working unchanged --
// tolerance values below 1 only ever distinguish "is the sign exactly
// +1" from "0 or -1", which is exactly the intended check.
double DelaunayTetrahedralization3D::inCircumsphere(const Tetrahedron& t, std::size_t pointIndex) const {
    const Vector3& p = points_[pointIndex];
    const Vector3& a = points_[t.vertices[0]];
    const Vector3& b = points_[t.vertices[1]];
    const Vector3& c = points_[t.vertices[2]];
    const Vector3& d = points_[t.vertices[3]];
    return static_cast<double>(inSphere3D(a, b, c, d, p));
}

namespace {

// Everything below is scratch state local to tetrahedralize()'s main loop
// -- never stored on the class, never touching the public Tetrahedron
// struct. It exists to answer, in O(1) amortized instead of O(current
// tetrahedron count), the one question findCavity()/findOtherFaceOwner()
// used to answer by scanning every tetrahedron: "which tetrahedron, if
// any, shares this face?"
constexpr std::size_t kNoNeighbor = std::numeric_limits<std::size_t>::max();

// One tetrahedron plus incremental adjacency. neighbor[li] mirrors
// tetrahedronFaces()'s own convention: the tetrahedron sharing the face
// opposite local vertex li, or kNoNeighbor if that face is currently
// exposed (owned by only one live tetrahedron so far). Dead tetrahedra are
// tombstoned (alive=false), never erased -- so no index anywhere ever
// needs correcting because some other element "moved".
struct WorkTet {
    DelaunayTetrahedralization3D::Tetrahedron t;
    std::array<std::size_t, 4> neighbor;
    bool alive;
};

// Faces owned by exactly one live WorkTet so far, keyed by sortedFace().
// A face with two owners is fully resolved (each side's `neighbor` already
// points at the other) and has no entry here. std::map, not
// unordered_map: avoids writing/testing a hash for array<size_t,3>, and
// this map's size tracks the perimeter of the currently-exposed region
// (small), not the point count, so O(log k) here is irrelevant next to
// the O(N) scans it replaces.
using FaceOwnerMap = std::map<std::array<std::size_t, 3>, std::pair<std::size_t, int>>;

// Registers tetIdx's local face li: links it to a previously-exposed owner
// on the other side if one exists (resolving the face -- it stops being
// an `owners` entry), or exposes it for the first time otherwise. Same
// "first touch exposes, second touch resolves" rule
// TetrahedralMesh::fromCells() already uses in one batch pass over a
// finished cell list -- mirrored here, applied incrementally instead.
void linkOrExposeFace(std::vector<WorkTet>& work, FaceOwnerMap& owners, std::size_t tetIdx, int li) {
    const Face f = tetrahedronFaces(work[tetIdx].t)[li];
    const auto key = sortedFace({f.a, f.b, f.c});
    const auto it = owners.find(key);
    if (it == owners.end()) {
        owners.emplace(key, std::make_pair(tetIdx, li));
        return;
    }
    const std::size_t otherIdx = it->second.first;
    const int otherLi = it->second.second;
    work[tetIdx].neighbor[static_cast<std::size_t>(li)] = otherIdx;
    work[otherIdx].neighbor[static_cast<std::size_t>(otherLi)] = tetIdx;
    owners.erase(it);
}

// Marks every tetrahedron in `bad` dead. Deliberately separate from the
// neighbor-fixup pass below: that pass must see the *final* liveness of
// both sides of a shared face, which is only guaranteed once every
// tetrahedron in this batch has already been marked -- doing both in one
// combined pass would let a face between two tetrahedra that are *both*
// dying in this same cavity replacement get wrongly re-exposed, depending
// on which of the two happened to be visited first.
void markDead(std::vector<WorkTet>& work, const std::vector<std::size_t>& bad) {
    for (std::size_t idx : bad) {
        work[idx].alive = false;
    }
}

// Re-exposes whichever faces of the now-dead tetrahedra in `bad` still
// have a *surviving* neighbor on the other side -- the inverse of
// linkOrExposeFace(). Must run after every element of `bad` has already
// had markDead() applied (see its own comment for why).
void reexposeSurvivingNeighbors(std::vector<WorkTet>& work, FaceOwnerMap& owners,
                                 const std::vector<std::size_t>& bad) {
    for (std::size_t deadIdx : bad) {
        for (int li = 0; li < 4; ++li) {
            const std::size_t nb = work[deadIdx].neighbor[static_cast<std::size_t>(li)];
            if (nb == kNoNeighbor) {
                const Face f = tetrahedronFaces(work[deadIdx].t)[li];
                owners.erase(sortedFace({f.a, f.b, f.c}));
                continue;
            }
            if (!work[nb].alive) {
                continue; // both sides of this internal face are dying together
            }
            for (int nbLi = 0; nbLi < 4; ++nbLi) {
                if (work[nb].neighbor[static_cast<std::size_t>(nbLi)] == deadIdx) {
                    work[nb].neighbor[static_cast<std::size_t>(nbLi)] = kNoNeighbor;
                    linkOrExposeFace(work, owners, nb, nbLi);
                    break;
                }
            }
        }
    }
}

// Point location: walks from `hint` toward point `p` across face
// adjacency, using isFaceVisible() (the same substitution test
// pointInTetrahedron() already uses) to decide which face, if any, p is
// on the far side of. Returns the (live) tetrahedron containing p, which
// is therefore guaranteed "bad" for Bowyer-Watson purposes: a
// tetrahedron's own interior always lies inside its circumsphere, since
// its four vertices lie exactly on that sphere and a ball is convex.
//
// Never revisits a tetrahedron within the same call -- tracked via
// `visitedStamp`/`generation` (a persistent array reused across every
// point, stamped rather than re-cleared, so this costs O(1) amortized to
// set up per call instead of O(current tetrahedron count)) -- which makes
// the walk self-bounding at `work.size()` steps regardless of geometry,
// closing off the classic failure mode of a naive "first violated face"
// visibility walk cycling forever. Returns kNoNeighbor only if every
// tetrahedron reachable from `hint` without revisiting is exhausted
// without finding one containing p: not expected for a point interior to
// the super-tetrahedron (every point tetrahedralize() inserts is, by
// construction), but the caller falls back to an exhaustive scan for that
// one point rather than trust the walk unconditionally.
std::size_t locateContaining(std::vector<WorkTet>& work, const std::vector<Vector3>& points,
                              std::vector<std::size_t>& visitedStamp, std::size_t generation,
                              std::size_t hint, std::size_t p) {
    std::size_t current = hint;
    while (true) {
        visitedStamp[current] = generation;
        bool anyViolated = false;
        std::size_t next = kNoNeighbor;
        for (int li = 0; li < 4; ++li) {
            if (!isFaceVisible(work[current].t, li, p, points)) {
                continue;
            }
            anyViolated = true;
            const std::size_t candidate = work[current].neighbor[static_cast<std::size_t>(li)];
            if (candidate != kNoNeighbor && visitedStamp[candidate] != generation) {
                next = candidate;
                break;
            }
        }
        if (!anyViolated) {
            return current; // p is inside (or on the boundary of) this tetrahedron
        }
        if (next == kNoNeighbor) {
            return kNoNeighbor; // every violated direction is a hull face or already visited
        }
        current = next;
    }
}

} // namespace

void DelaunayTetrahedralization3D::tetrahedralize() {
    const std::size_t originalPointCount = points_.size();
    tetrahedra_.clear();
    if (originalPointCount == 0) {
        return;
    }

    // Bounding box and super-tetrahedron sizing: identical to
    // tetrahedralizeReference() -- see its own extensive comment for why R
    // must be this large. Duplicated rather than shared because the two
    // functions otherwise share no state (one builds `tetrahedra_`
    // directly, this one builds a local `work` list), and the block is
    // small enough that a shared helper would need to return six values
    // for little benefit.
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    double minZ = std::numeric_limits<double>::max();
    double maxZ = std::numeric_limits<double>::lowest();
    for (std::size_t i = 0; i < originalPointCount; ++i) {
        minX = std::min(minX, points_[i].x);
        maxX = std::max(maxX, points_[i].x);
        minY = std::min(minY, points_[i].y);
        maxY = std::max(maxY, points_[i].y);
        minZ = std::min(minZ, points_[i].z);
        maxZ = std::max(maxZ, points_[i].z);
    }
    const double dx = maxX - minX;
    const double dy = maxY - minY;
    const double dz = maxZ - minZ;
    const double deltaMax = std::max({dx, dy, dz}) + 1.0;
    const double midX = (minX + maxX) / 2.0;
    const double midY = (minY + maxY) / 2.0;
    const double midZ = (minZ + maxZ) / 2.0;
    const double r = 1.0e8 * deltaMax;
    const std::size_t superA = points_.size();
    points_.emplace_back(midX + r, midY + r, midZ + r);
    const std::size_t superB = points_.size();
    points_.emplace_back(midX + r, midY - r, midZ - r);
    const std::size_t superC = points_.size();
    points_.emplace_back(midX - r, midY + r, midZ - r);
    const std::size_t superD = points_.size();
    points_.emplace_back(midX - r, midY - r, midZ + r);

    std::vector<WorkTet> work;
    work.reserve(originalPointCount * 4 + 8);
    std::vector<std::size_t> visitedStamp; // point-location walk's own visited marker
    visitedStamp.reserve(originalPointCount * 4 + 8);
    std::vector<std::size_t> badStamp; // cavity BFS's own "in this point's bad set" marker
    badStamp.reserve(originalPointCount * 4 + 8);
    FaceOwnerMap owners;

    work.push_back(WorkTet{makeTetrahedron(superA, superB, superC, superD, points_),
                            {kNoNeighbor, kNoNeighbor, kNoNeighbor, kNoNeighbor}, true});
    visitedStamp.push_back(0);
    badStamp.push_back(0);
    for (int li = 0; li < 4; ++li) {
        linkOrExposeFace(work, owners, 0, li);
    }

    std::size_t hint = 0;
    std::vector<std::size_t> badIndices;
    std::vector<std::array<std::size_t, 3>> boundary;
    std::vector<std::size_t> queue;

    for (std::size_t p = 0; p < originalPointCount; ++p) {
        const std::size_t generation = p + 1; // 0 is reserved as "never visited/never bad"

        std::size_t seedTet = locateContaining(work, points_, visitedStamp, generation, hint, p);
        if (seedTet == kNoNeighbor) {
            // Not expected for a point interior to the super-tetrahedron
            // (see locateContaining()'s own comment) -- correctness must
            // never depend on the walk succeeding, only speed does. Falls
            // back to an exhaustive scan of live tetrahedra for this one
            // point; some live tetrahedron's circumsphere is guaranteed to
            // contain p, by the same interior-to-the-super-tetrahedron
            // argument tetrahedralizeReference() itself relies on.
            for (std::size_t t = 0; t < work.size(); ++t) {
                if (work[t].alive && inCircumsphere(work[t].t, p) > 0.0) {
                    seedTet = t;
                    break;
                }
            }
        }

        badIndices.clear();
        boundary.clear();
        queue.clear();
        queue.push_back(seedTet);
        badIndices.push_back(seedTet);
        badStamp[seedTet] = generation;

        std::size_t qi = 0;
        while (qi < queue.size()) {
            const std::size_t cur = queue[qi++];
            for (int li = 0; li < 4; ++li) {
                const std::size_t nb = work[cur].neighbor[static_cast<std::size_t>(li)];
                const Face f = tetrahedronFaces(work[cur].t)[li];
                if (nb == kNoNeighbor) {
                    // Hull face of the current live region: unconditionally
                    // boundary, same semantics findCavity() documents for
                    // considerHullVisibility=false -- always correct here
                    // since every point tetrahedralize() inserts is
                    // strictly interior to the enclosing super-tetrahedron.
                    boundary.push_back({f.a, f.b, f.c});
                    continue;
                }
                if (badStamp[nb] == generation) {
                    continue; // interior to the cavity, already counted from the other side
                }
                if (inCircumsphere(work[nb].t, p) > 0.0) {
                    badStamp[nb] = generation;
                    badIndices.push_back(nb);
                    queue.push_back(nb);
                } else {
                    boundary.push_back({f.a, f.b, f.c});
                }
            }
        }

        markDead(work, badIndices);
        reexposeSurvivingNeighbors(work, owners, badIndices);

        for (const auto& f : boundary) {
            WorkTet nt;
            nt.t = makeTetrahedron(f[0], f[1], f[2], p, points_);
            nt.neighbor = {kNoNeighbor, kNoNeighbor, kNoNeighbor, kNoNeighbor};
            nt.alive = true;
            work.push_back(nt);
            visitedStamp.push_back(0);
            badStamp.push_back(0);
            const std::size_t newIdx = work.size() - 1;
            for (int li = 0; li < 4; ++li) {
                linkOrExposeFace(work, owners, newIdx, li);
            }
            hint = newIdx;
        }
    }

    // Compaction: the same "drop anything touching a super vertex" filter
    // tetrahedralizeReference() applies at the end, with one extra
    // criterion -- liveness -- since dead tetrahedra were tombstoned
    // rather than physically removed during the loop above.
    std::vector<Tetrahedron> result;
    result.reserve(work.size());
    for (const WorkTet& wt : work) {
        if (!wt.alive) {
            continue;
        }
        if (tetrahedronHasVertex(wt.t, superA) || tetrahedronHasVertex(wt.t, superB) ||
            tetrahedronHasVertex(wt.t, superC) || tetrahedronHasVertex(wt.t, superD)) {
            continue;
        }
        result.push_back(wt.t);
    }
    tetrahedra_ = std::move(result);

    points_.resize(originalPointCount);
}

void DelaunayTetrahedralization3D::tetrahedralizeReference() {
    const std::size_t originalPointCount = points_.size();
    tetrahedra_.clear();
    if (originalPointCount == 0) {
        return;
    }

    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    double minZ = std::numeric_limits<double>::max();
    double maxZ = std::numeric_limits<double>::lowest();
    for (std::size_t i = 0; i < originalPointCount; ++i) {
        minX = std::min(minX, points_[i].x);
        maxX = std::max(maxX, points_[i].x);
        minY = std::min(minY, points_[i].y);
        maxY = std::max(maxY, points_[i].y);
        minZ = std::min(minZ, points_[i].z);
        maxZ = std::max(maxZ, points_[i].z);
    }
    const double dx = maxX - minX;
    const double dy = maxY - minY;
    const double dz = maxZ - minZ;
    const double deltaMax = std::max({dx, dy, dz}) + 1.0; // +1 guards a degenerate single-point case
    const double midX = (minX + maxX) / 2.0;
    const double midY = (minY + maxY) / 2.0;
    const double midZ = (minZ + maxZ) / 2.0;

    // Super tetrahedron: alternating corners of a cube of half-side R,
    // centered on the input's bounding box -- a regular tetrahedron whose
    // inradius (~0.577*R) comfortably exceeds the input's own extent from
    // the center. Its vertices are appended after the real points so their
    // indices are known, and every tetrahedron touching them is discarded
    // at the end.
    //
    // **R must be far larger than merely "big enough to enclose every
    // input point"** -- that weaker condition (satisfied from R ≈
    // deltaMax onward, and the reason an earlier version used only
    // 20*deltaMax) is genuinely not sufficient, and getting it wrong
    // caused a real bug: the final tetrahedralization occasionally came
    // out missing a thin, valid, hull-adjacent tetrahedron entirely.
    //
    // Why: the super-tetrahedron's four vertices are a *finite* stand-in
    // for "points at infinity". A tetrahedron built on some real points
    // plus a super vertex therefore gets a finite circumsphere where the
    // true (conceptual) one is a half-space. If the super vertices are
    // only moderately far away, that finite circumsphere can be small
    // enough to wrongly test as containing a nearby real point, deleting
    // a tetrahedron that should have survived -- and a sliver near the
    // convex hull is exactly where that error surfaces. Pushing the super
    // vertices much farther out makes those circumspheres approximate the
    // true half-spaces closely enough for the in-sphere test to answer
    // correctly.
    //
    // Confirmed by direct measurement (exact rational arithmetic, so
    // floating-point precision plays no part in these numbers): on the
    // 7-point set that first exposed this, R = 20*deltaMax produced 6
    // tetrahedra and a volume 0.00096 short of the true convex hull
    // volume, while R >= 200*deltaMax produced the correct 7 matching
    // qhull exactly. A later randomized sweep found a 12-point set still
    // failing at 2000*deltaMax (volume 6.1e-05 short) and correct from
    // 20000*deltaMax up.
    //
    // The scale below therefore carries a large safety margin over the
    // largest failure actually observed, rather than sitting just past
    // it. This is cheap: enlarging the super-tetrahedron costs nothing at
    // runtime (it is four points, discarded at the end), and the
    // coordinates stay far inside double's exponent range for any
    // realistic input, so there is no precision downside to being
    // generous here.
    //
    // Strictly speaking this is a *mitigation* rather than a proof: the
    // fully principled alternative is to treat the super vertices
    // symbolically as true points at infinity (an infinite cell's
    // "circumsphere" is then the half-space outside its hull face, so the
    // in-sphere test against it reduces exactly to an orientation test),
    // which needs no scale constant at all and is what production
    // implementations do.
    //
    // That alternative was prototyped and validated separately (80/80 on
    // adversarial inputs: near-coplanar slabs, co-spherical points, and
    // extreme scale disparity), then compared head-to-head against the
    // code below. **They agree exactly** -- identical tetrahedra and
    // bit-for-bit identical exact volumes on every case tried, including
    // clusters spanning 1e-3 with an outlier at 1e6. So the scale below
    // is not merely "large enough so far": it is empirically
    // indistinguishable from the scale-free formulation on everything
    // measurable. The symbolic rewrite was therefore deliberately not
    // adopted -- it would rewrite this class's core loop (recently and
    // painfully debugged) to fix no observed failure. Its one remaining
    // attraction is that it would also subsume findCavity()'s
    // considerHullVisibility special case, which is worth revisiting if
    // that path ever needs to grow.
    const double r = 1.0e8 * deltaMax;
    const std::size_t superA = points_.size();
    points_.emplace_back(midX + r, midY + r, midZ + r);
    const std::size_t superB = points_.size();
    points_.emplace_back(midX + r, midY - r, midZ - r);
    const std::size_t superC = points_.size();
    points_.emplace_back(midX - r, midY + r, midZ - r);
    const std::size_t superD = points_.size();
    points_.emplace_back(midX - r, midY - r, midZ + r);

    tetrahedra_.push_back(makeTetrahedron(superA, superB, superC, superD, points_));

    for (std::size_t p = 0; p < originalPointCount; ++p) {
        std::vector<std::size_t> badIndices;
        std::vector<std::array<std::size_t, 3>> boundary;
        findCavity(p, /*considerHullVisibility=*/false, badIndices, boundary);

        std::vector<Tetrahedron> kept;
        kept.reserve(tetrahedra_.size());
        for (std::size_t t = 0; t < tetrahedra_.size(); ++t) {
            if (std::find(badIndices.begin(), badIndices.end(), t) == badIndices.end()) {
                kept.push_back(tetrahedra_[t]);
            }
        }
        for (const auto& f : boundary) {
            kept.push_back(makeTetrahedron(f[0], f[1], f[2], p, points_));
        }
        tetrahedra_ = std::move(kept);
    }

    std::vector<Tetrahedron> withoutSuper;
    withoutSuper.reserve(tetrahedra_.size());
    for (const Tetrahedron& t : tetrahedra_) {
        if (!tetrahedronHasVertex(t, superA) && !tetrahedronHasVertex(t, superB) &&
            !tetrahedronHasVertex(t, superC) && !tetrahedronHasVertex(t, superD)) {
            withoutSuper.push_back(t);
        }
    }
    tetrahedra_ = std::move(withoutSuper);

    points_.resize(originalPointCount); // drop the super tetrahedron's vertices
}

bool DelaunayTetrahedralization3D::insertSteinerPoint(double x, double y, double z) {
    if (tetrahedra_.empty()) {
        return false; // nothing to insert into
    }

    const std::size_t p = points_.size();
    points_.emplace_back(x, y, z);

    std::vector<std::size_t> badIndices;
    std::vector<std::array<std::size_t, 3>> boundary;
    findCavity(p, /*considerHullVisibility=*/true, badIndices, boundary);

    // findCavity() finds the correct conflict region whether p is interior
    // or exterior to the current hull (see its own comment). For a purely
    // exterior p that doesn't invalidate any existing tetrahedron's own
    // circumsphere, badIndices legitimately comes back *empty* -- no real
    // tetrahedron needs removing, only new cap tetrahedra get added over
    // the visible hull faces (exactly like the earlier hand-derived
    // example this function's comment references). So `boundary` being
    // non-empty, not `badIndices`, is the real "is there anything to do"
    // signal -- for any point distinct from the existing set, boundary
    // should never come back empty (a convex hull always has at least one
    // face visible from any exterior point, and the Delaunay tetrahedra
    // already cover every interior point).
    //
    // Note the boundary faces are *not* expected to form a closed shell by
    // themselves when p is exterior: the discarded visible hull face(s)
    // leave a genuine "opening" that p's new tetrahedra fill, unlike the
    // interior case where the boundary is a closed cavity wall. So the
    // only defensive guard here is against a genuinely degenerate p (e.g.
    // exactly coincident with an existing point, or an exact coplanarity/
    // co-sphericity tie): reject if any resulting new tetrahedron would be
    // degenerate (near-zero volume), rather than build an invalid
    // tetrahedralization.
    if (boundary.empty()) {
        points_.pop_back();
        return false;
    }

    // Exact zero (not a tolerance-based "close to zero") now that this is
    // orientation3D's exact sign, not a plain double computation -- a
    // genuinely coplanar/degenerate case gives exactly 0, never a small
    // nonzero value to have to guess a threshold for.
    for (const auto& f : boundary) {
        const Vector3& a = points_[f[0]];
        const Vector3& b = points_[f[1]];
        const Vector3& c = points_[f[2]];
        const Vector3& d = points_[p];
        if (orientation3D(a, b, c, d) == 0) {
            points_.pop_back();
            return false;
        }
    }

    std::vector<Tetrahedron> kept;
    kept.reserve(tetrahedra_.size());
    for (std::size_t t = 0; t < tetrahedra_.size(); ++t) {
        if (std::find(badIndices.begin(), badIndices.end(), t) == badIndices.end()) {
            kept.push_back(tetrahedra_[t]);
        }
    }
    for (const auto& f : boundary) {
        kept.push_back(makeTetrahedron(f[0], f[1], f[2], p, points_));
    }
    tetrahedra_ = std::move(kept);
    return true;
}

bool DelaunayTetrahedralization3D::hasFace(const std::array<std::size_t, 3>& face) const {
    const auto key = sortedFace(face);
    for (const Tetrahedron& t : tetrahedra_) {
        for (const Face& f : tetrahedronFaces(t)) {
            if (sortedFace({f.a, f.b, f.c}) == key) {
                return true;
            }
        }
    }
    return false;
}

std::vector<std::array<std::size_t, 3>> DelaunayTetrahedralization3D::missingFacets(
    const std::vector<std::array<std::size_t, 3>>& facets) const {
    std::set<std::array<std::size_t, 3>> present;
    for (const Tetrahedron& t : tetrahedra_) {
        for (const Face& f : tetrahedronFaces(t)) {
            present.insert(sortedFace({f.a, f.b, f.c}));
        }
    }
    std::vector<std::array<std::size_t, 3>> missing;
    for (const auto& f : facets) {
        if (present.find(sortedFace(f)) == present.end()) {
            missing.push_back(f);
        }
    }
    return missing;
}

void DelaunayTetrahedralization3D::recoverFacetRecursive(const std::array<std::size_t, 3>& facet,
                                                           int roundsLeft,
                                                           std::vector<std::array<std::size_t, 3>>& recovered,
                                                           std::vector<std::array<std::size_t, 3>>& unrecovered) {
    if (hasFace(facet)) {
        recovered.push_back(facet);
        return;
    }
    if (roundsLeft <= 0) {
        unrecovered.push_back(facet);
        return;
    }
    const Vector3 centroid = (points_[facet[0]] + points_[facet[1]] + points_[facet[2]]) / 3.0;
    if (!insertSteinerPoint(centroid.x, centroid.y, centroid.z)) {
        unrecovered.push_back(facet);
        return;
    }
    const std::size_t c = points_.size() - 1;
    recoverFacetRecursive({facet[0], facet[1], c}, roundsLeft - 1, recovered, unrecovered);
    recoverFacetRecursive({facet[1], facet[2], c}, roundsLeft - 1, recovered, unrecovered);
    recoverFacetRecursive({facet[2], facet[0], c}, roundsLeft - 1, recovered, unrecovered);
}

bool DelaunayTetrahedralization3D::tryFlipCoplanarQuadDiagonal(const std::array<std::size_t, 3>& f1,
                                                                 const std::array<std::size_t, 3>& f2) {
    // Identify the shared edge {sharedA, sharedB} = f1 ∩ f2 (the diagonal
    // this flip would introduce) and each triangle's own vertex not in the
    // other (q from f1, s from f2 -- together with sharedA/sharedB, the
    // quad's other diagonal). Anything but exactly 2 shared vertices means
    // f1/f2 do not describe two triangles of one shared quadrilateral.
    std::array<std::size_t, 2> shared{};
    std::size_t sharedCount = 0;
    std::size_t q = 0;
    bool qSet = false;
    for (std::size_t v : f1) {
        const bool inF2 = (v == f2[0] || v == f2[1] || v == f2[2]);
        if (inF2) {
            if (sharedCount >= 2) {
                return false; // f1 == f2 as a vertex set
            }
            shared[sharedCount++] = v;
        } else {
            q = v;
            qSet = true;
        }
    }
    if (sharedCount != 2 || !qSet) {
        return false;
    }
    std::size_t s = 0;
    bool sSet = false;
    for (std::size_t v : f2) {
        const bool inF1 = (v == f1[0] || v == f1[1] || v == f1[2]);
        if (!inF1) {
            if (sSet) {
                return false; // more than one vertex of f2 absent from f1
            }
            s = v;
            sSet = true;
        }
    }
    if (!sSet || s == q) {
        return false;
    }
    const std::size_t sharedA = shared[0];
    const std::size_t sharedB = shared[1];

    // Coplanarity: sharedA, sharedB, q and s must lie in one plane --
    // otherwise f1/f2 do not describe a single planar quadrilateral and no
    // diagonal flip relates them.
    if (orientation3D(points_[sharedA], points_[sharedB], points_[q], points_[s]) != 0) {
        return false;
    }

    // The candidate quad, in cyclic order, is sharedA-q-sharedB-s (edge
    // sharedA-q and q-sharedB from f1, sharedB-s and s-sharedA from f2).
    // {sharedA, sharedB} is the diagonal this flip introduces; {q, s} is
    // the one the tetrahedralization currently uses instead, which must
    // presently split the quad into exactly two tetrahedra sharing one
    // common fourth vertex -- the simple "bipyramid" case this flip
    // handles. A different local structure (the two triangles owned by
    // tetrahedra with different fourth vertices, or missing outright) is
    // left to the caller's Steiner fallback.
    const FaceOwner ownerQ = findFaceOwnerWithApex(tetrahedra_, {sharedA, q, s});
    const FaceOwner ownerS = findFaceOwnerWithApex(tetrahedra_, {sharedB, q, s});
    if (ownerQ.tetIndex < 0 || ownerS.tetIndex < 0 || ownerQ.apex != ownerS.apex) {
        return false;
    }
    const std::size_t apex = ownerQ.apex;

    // Convexity of quad sharedA-q-sharedB-s, seen consistently from `apex`:
    // all four consecutive-triple orientations must agree in sign, the
    // same turn-direction test a 2D convex-polygon check uses, with `apex`
    // fixing which side of the shared plane counts as positive. A
    // non-convex quad would make this diagonal an invalid split (the new
    // tetrahedra would overlap or leave a gap instead of exactly
    // partitioning the same volume the old two did) -- checked explicitly
    // here rather than trusted, since makeTetrahedron() below would
    // "succeed" at building a positively-oriented tetrahedron regardless of
    // whether the split is geometrically sound.
    const int t1 = orientation3D(points_[sharedA], points_[q], points_[sharedB], points_[apex]);
    const int t2 = orientation3D(points_[q], points_[sharedB], points_[s], points_[apex]);
    const int t3 = orientation3D(points_[sharedB], points_[s], points_[sharedA], points_[apex]);
    const int t4 = orientation3D(points_[s], points_[sharedA], points_[q], points_[apex]);
    if (t1 == 0 || t2 == 0 || t3 == 0 || t4 == 0 || !(t1 == t2 && t2 == t3 && t3 == t4)) {
        return false;
    }

    const auto t1Index = static_cast<std::size_t>(ownerQ.tetIndex);
    const auto t2Index = static_cast<std::size_t>(ownerS.tetIndex);
    std::vector<Tetrahedron> kept;
    kept.reserve(tetrahedra_.size());
    for (std::size_t t = 0; t < tetrahedra_.size(); ++t) {
        if (t != t1Index && t != t2Index) {
            kept.push_back(tetrahedra_[t]);
        }
    }
    // {sharedA, q, sharedB, apex} exposes face f1 = {sharedA, q, sharedB};
    // {sharedA, sharedB, s, apex} exposes face f2 = {sharedA, sharedB, s}.
    // Together they occupy exactly the same volume the two removed
    // tetrahedra did (same base quad, same apex, only the internal
    // diagonal differs), confirmed convex above.
    kept.push_back(makeTetrahedron(sharedA, q, sharedB, apex, points_));
    kept.push_back(makeTetrahedron(sharedA, sharedB, s, apex, points_));
    tetrahedra_ = std::move(kept);
    return true;
}

DelaunayTetrahedralization3D::FacetRecoveryResult DelaunayTetrahedralization3D::recoverFacets(
    const std::vector<std::array<std::size_t, 3>>& facets, int maxRounds) {
    // Phase 1: combinatorial recovery, tried before any new point --
    // repeatedly scan the still-missing facets for a flippable coplanar-quad
    // pair (see tryFlipCoplanarQuadDiagonal()) until a full pass finds none.
    // Bounded by construction: each successful flip strictly reduces the
    // missing count by 2, so this terminates in at most facets.size()/2
    // successful iterations regardless of the outer restart-on-success
    // pattern below.
    std::vector<std::array<std::size_t, 3>> stillMissing = missingFacets(facets);
    std::vector<bool> flipped(stillMissing.size(), false);
    bool progress = true;
    while (progress) {
        progress = false;
        for (std::size_t i = 0; i < stillMissing.size() && !progress; ++i) {
            if (flipped[i]) {
                continue;
            }
            for (std::size_t j = i + 1; j < stillMissing.size(); ++j) {
                if (flipped[j]) {
                    continue;
                }
                if (tryFlipCoplanarQuadDiagonal(stillMissing[i], stillMissing[j])) {
                    flipped[i] = true;
                    flipped[j] = true;
                    progress = true;
                    break;
                }
            }
        }
    }

    // Phase 2: every facet in the original request, in its original order.
    // hasFace() (recoverFacetRecursive()'s first check) already reflects
    // phase 1's flips, so a facet fixed above costs nothing extra here; an
    // untouched facet falls through to the centroid-Steiner heuristic
    // exactly as before phase 1 existed.
    FacetRecoveryResult result;
    for (const auto& f : facets) {
        recoverFacetRecursive(f, maxRounds, result.recoveredFacets, result.unrecovered);
    }
    return result;
}

std::size_t DelaunayTetrahedralization3D::removeRegion(const Vector3& seed,
                                                         const std::vector<std::array<std::size_t, 3>>& walls) {
    long long seedTet = -1;
    for (std::size_t t = 0; t < tetrahedra_.size(); ++t) {
        if (pointInTetrahedron(tetrahedra_[t], seed, points_)) {
            seedTet = static_cast<long long>(t);
            break;
        }
    }
    if (seedTet == -1) {
        return 0;
    }

    std::set<std::array<std::size_t, 3>> wallKeys;
    for (const auto& w : walls) {
        wallKeys.insert(sortedFace(w));
    }

    std::vector<bool> visited(tetrahedra_.size(), false);
    std::vector<std::size_t> stack = {static_cast<std::size_t>(seedTet)};
    visited[static_cast<std::size_t>(seedTet)] = true;
    std::vector<std::size_t> toRemove;

    while (!stack.empty()) {
        const std::size_t current = stack.back();
        stack.pop_back();
        toRemove.push_back(current);
        for (const Face& f : tetrahedronFaces(tetrahedra_[current])) {
            if (wallKeys.find(sortedFace({f.a, f.b, f.c})) != wallKeys.end()) {
                continue; // impassable recovered boundary
            }
            const long long other = findOtherFaceOwner(tetrahedra_, f, current);
            if (other == -1 || visited[static_cast<std::size_t>(other)]) {
                continue;
            }
            visited[static_cast<std::size_t>(other)] = true;
            stack.push_back(static_cast<std::size_t>(other));
        }
    }

    std::sort(toRemove.begin(), toRemove.end());
    std::vector<Tetrahedron> kept;
    kept.reserve(tetrahedra_.size() - toRemove.size());
    std::size_t ri = 0;
    for (std::size_t t = 0; t < tetrahedra_.size(); ++t) {
        if (ri < toRemove.size() && toRemove[ri] == t) {
            ++ri;
            continue;
        }
        kept.push_back(tetrahedra_[t]);
    }
    tetrahedra_ = std::move(kept);
    return toRemove.size();
}

bool DelaunayTetrahedralization3D::satisfiesDelaunayProperty(double tolerance) const {
    for (const Tetrahedron& t : tetrahedra_) {
        for (std::size_t p = 0; p < points_.size(); ++p) {
            if (tetrahedronHasVertex(t, p)) {
                continue;
            }
            if (inCircumsphere(t, p) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

} // namespace aether::mesh
