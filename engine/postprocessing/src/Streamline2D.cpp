#include "aether/postprocessing/Streamline2D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::postprocessing {

using aether::core::Vector3;

Streamline2D::Streamline2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY,
                            const std::vector<double>& u, const std::vector<double>& v, bool periodic)
    : nx_(nx), ny_(ny), lengthX_(lengthX), lengthY_(lengthY), dx_(lengthX / static_cast<double>(nx)),
      dy_(lengthY / static_cast<double>(ny)), periodic_(periodic), u_(u), v_(v) {}

namespace {

double sampleField(const std::vector<double>& field, std::size_t nx, std::size_t ny, long long i,
                    long long j, bool periodic) {
    if (periodic) {
        i = ((i % static_cast<long long>(nx)) + static_cast<long long>(nx)) % static_cast<long long>(nx);
        j = ((j % static_cast<long long>(ny)) + static_cast<long long>(ny)) % static_cast<long long>(ny);
    } else {
        i = std::clamp(i, 0LL, static_cast<long long>(nx) - 1);
        j = std::clamp(j, 0LL, static_cast<long long>(ny) - 1);
    }
    return field[static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * nx];
}

} // namespace

Vector3 Streamline2D::velocityAt(double x, double y) const {
    const double fx = x / dx_ - 0.5;
    const double fy = y / dy_ - 0.5;
    const auto i0 = static_cast<long long>(std::floor(fx));
    const auto j0 = static_cast<long long>(std::floor(fy));
    const double tx = fx - static_cast<double>(i0);
    const double ty = fy - static_cast<double>(j0);
    const long long i1 = i0 + 1;
    const long long j1 = j0 + 1;

    const double u00 = sampleField(u_, nx_, ny_, i0, j0, periodic_);
    const double u10 = sampleField(u_, nx_, ny_, i1, j0, periodic_);
    const double u01 = sampleField(u_, nx_, ny_, i0, j1, periodic_);
    const double u11 = sampleField(u_, nx_, ny_, i1, j1, periodic_);
    const double uInterp =
        (1.0 - tx) * (1.0 - ty) * u00 + tx * (1.0 - ty) * u10 + (1.0 - tx) * ty * u01 + tx * ty * u11;

    const double v00 = sampleField(v_, nx_, ny_, i0, j0, periodic_);
    const double v10 = sampleField(v_, nx_, ny_, i1, j0, periodic_);
    const double v01 = sampleField(v_, nx_, ny_, i0, j1, periodic_);
    const double v11 = sampleField(v_, nx_, ny_, i1, j1, periodic_);
    const double vInterp =
        (1.0 - tx) * (1.0 - ty) * v00 + tx * (1.0 - ty) * v10 + (1.0 - tx) * ty * v01 + tx * ty * v11;

    return Vector3(uInterp, vInterp, 0.0);
}

std::vector<Vector3> Streamline2D::trace(double x0, double y0, double stepSize, std::size_t maxSteps) const {
    std::vector<Vector3> path;
    path.reserve(maxSteps + 1);
    double x = x0;
    double y = y0;
    path.emplace_back(x, y, 0.0);

    for (std::size_t step = 0; step < maxSteps; ++step) {
        const Vector3 k1 = velocityAt(x, y);
        if (k1.normSquared() < 1e-24) {
            break; // stagnation point
        }
        const Vector3 k2 = velocityAt(x + 0.5 * stepSize * k1.x, y + 0.5 * stepSize * k1.y);
        const Vector3 k3 = velocityAt(x + 0.5 * stepSize * k2.x, y + 0.5 * stepSize * k2.y);
        const Vector3 k4 = velocityAt(x + stepSize * k3.x, y + stepSize * k3.y);

        double newX = x + (stepSize / 6.0) * (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x);
        double newY = y + (stepSize / 6.0) * (k1.y + 2.0 * k2.y + 2.0 * k3.y + k4.y);

        if (periodic_) {
            newX = std::fmod(newX, lengthX_);
            if (newX < 0.0) {
                newX += lengthX_;
            }
            newY = std::fmod(newY, lengthY_);
            if (newY < 0.0) {
                newY += lengthY_;
            }
        } else if (newX < 0.0 || newX > lengthX_ || newY < 0.0 || newY > lengthY_) {
            break; // left the domain
        }

        x = newX;
        y = newY;
        path.emplace_back(x, y, 0.0);
    }
    return path;
}

} // namespace aether::postprocessing
