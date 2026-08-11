#include "aether/core/Matrix4x4.hpp"
#include "aether/core/Mesh.hpp"
#include "aether/core/ScalarField.hpp"
#include "aether/core/Tensor3x3.hpp"
#include "aether/core/Vector3.hpp"
#include "aether/core/VectorField.hpp"
#include "aether/testing/Check.hpp"

#include <cmath>
#include <cstdio>

using namespace aether::core;

namespace {

bool nearlyEqual(double a, double b, double tolerance = 1e-12) {
    return std::fabs(a - b) <= tolerance;
}

void testVector3() {
    const Vector3 a(1.0, 2.0, 3.0);
    const Vector3 b(4.0, 5.0, 6.0);

    AETHER_CHECK(nearlyEqual((a + b).x, 5.0));
    AETHER_CHECK(nearlyEqual(a.dot(b), 32.0));

    const Vector3 cross = a.cross(b);
    AETHER_CHECK(nearlyEqual(cross.x, -3.0));
    AETHER_CHECK(nearlyEqual(cross.y, 6.0));
    AETHER_CHECK(nearlyEqual(cross.z, -3.0));

    const Vector3 unit(3.0, 0.0, 4.0);
    AETHER_CHECK(nearlyEqual(unit.norm(), 5.0));
    AETHER_CHECK(nearlyEqual(unit.normalized().norm(), 1.0));
}

void testTensor3x3() {
    const Tensor3x3 identity = Tensor3x3::identity();
    const Vector3 v(1.0, 2.0, 3.0);
    const Vector3 result = identity * v;
    AETHER_CHECK(result == v);

    const Tensor3x3 a(1, 2, 3, 4, 5, 6, 7, 8, 9);
    AETHER_CHECK(nearlyEqual(a.trace(), 15.0));
    AETHER_CHECK(nearlyEqual(a.transposed()(0, 1), a(1, 0)));
}

void testMesh() {
    Mesh mesh;
    const std::size_t v0 = mesh.addVertex({0.0, 0.0, 0.0});
    const std::size_t v1 = mesh.addVertex({2.0, 0.0, 0.0});
    const std::size_t v2 = mesh.addVertex({0.0, 2.0, 0.0});
    mesh.addCell({v0, v1, v2});

    AETHER_CHECK(mesh.vertexCount() == 3);
    AETHER_CHECK(mesh.cellCount() == 1);

    const Vector3 centroid = mesh.cellCentroid(0);
    AETHER_CHECK(nearlyEqual(centroid.x, 2.0 / 3.0));
    AETHER_CHECK(nearlyEqual(centroid.y, 2.0 / 3.0));
}

// Self-derived, not a recalled reference matrix: for eye=(0,0,5),
// center=(0,0,0), up=(0,1,0), the forward axis f=(0,0,-1), side s=(1,0,0),
// up'=(0,1,0) -- so the origin (the look-at target) must land at camera-
// space (0,0,-5): the standard right-handed-view-space convention (camera
// looks down -z), at exactly the eye-to-target distance.
void testMatrix4x4LookAtPlacesTargetAtNegativeDistance() {
    const Matrix4x4 view = Matrix4x4::lookAt(Vector3(0.0, 0.0, 5.0), Vector3(0.0, 0.0, 0.0),
                                              Vector3(0.0, 1.0, 0.0));
    const Vector3 originInView = transformPoint(view, Vector3(0.0, 0.0, 0.0));
    AETHER_CHECK(nearlyEqual(originInView.x, 0.0, 1e-6));
    AETHER_CHECK(nearlyEqual(originInView.y, 0.0, 1e-6));
    AETHER_CHECK(nearlyEqual(originInView.z, -5.0, 1e-6));
}

// Self-derived: ortho(0,2,0,1,-1,1) maps its domain's exact center
// (1, 0.5, 0) to the NDC origin (0,0,0) -- true for any orthographic
// projection by construction (it is an affine map centered on the box's
// midpoint), not something that needs an external reference value.
void testMatrix4x4OrthoMapsBoxCenterToOrigin() {
    const Matrix4x4 proj = Matrix4x4::ortho(0.0, 2.0, 0.0, 1.0, -1.0, 1.0);
    const Vector3 center = transformPoint(proj, Vector3(1.0, 0.5, 0.0));
    AETHER_CHECK(nearlyEqual(center.x, 0.0, 1e-6));
    AETHER_CHECK(nearlyEqual(center.y, 0.0, 1e-6));
    AETHER_CHECK(nearlyEqual(center.z, 0.0, 1e-6));
}

void testMatrix4x4IdentityIsMultiplicativeIdentity() {
    const Matrix4x4 view = Matrix4x4::lookAt(Vector3(1.0, 2.0, 3.0), Vector3(0.0, 0.0, 0.0),
                                              Vector3(0.0, 1.0, 0.0));
    const Matrix4x4 product = Matrix4x4::identity() * view;
    for (std::size_t i = 0; i < 16; ++i) {
        AETHER_CHECK(nearlyEqual(product.m[i], view.m[i], 1e-9));
    }
}

void testFields() {
    Mesh mesh;
    const std::size_t v0 = mesh.addVertex({0.0, 0.0, 0.0});
    const std::size_t v1 = mesh.addVertex({1.0, 0.0, 0.0});
    const std::size_t v2 = mesh.addVertex({0.0, 1.0, 0.0});
    mesh.addCell({v0, v1, v2});
    mesh.addCell({v0, v1, v2});

    ScalarField pressure(mesh, 1.0);
    AETHER_CHECK(pressure.size() == 2);
    pressure[0] = 5.0;
    const ScalarField doubled = pressure * 2.0;
    AETHER_CHECK(nearlyEqual(doubled[0], 10.0));
    AETHER_CHECK(nearlyEqual(doubled[1], 2.0));

    VectorField velocity(mesh, Vector3(1.0, 0.0, 0.0));
    const VectorField sum = velocity + velocity;
    AETHER_CHECK(nearlyEqual(sum[0].x, 2.0));
}

} // namespace

int main() {
    testVector3();
    testTensor3x3();
    testMesh();
    testFields();
    testMatrix4x4LookAtPlacesTargetAtNegativeDistance();
    testMatrix4x4OrthoMapsBoxCenterToOrigin();
    testMatrix4x4IdentityIsMultiplicativeIdentity();
    std::puts("aether_core_tests: all tests passed");
    return 0;
}
