#include "aether/solver/CompressibleEulerSolver1D.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace aether::solver {

CompressibleEulerSolver1D::CompressibleEulerSolver1D(std::size_t n, double length, double gamma,
                                                       BoundaryCondition boundary)
    : n_(n), length_(length), dx_(length / static_cast<double>(n)), gamma_(gamma), boundary_(boundary),
      rho_(n, 1.0), rhoU_(n, 0.0), E_(n, 1.0) {}

void CompressibleEulerSolver1D::initialize(
    const std::function<void(double, double&, double&, double&)>& profile) {
    for (std::size_t i = 0; i < n_; ++i) {
        double rho = 0.0;
        double u = 0.0;
        double p = 0.0;
        profile(cellCenter(i), rho, u, p);
        rho_[i] = rho;
        rhoU_[i] = rho * u;
        E_[i] = p / (gamma_ - 1.0) + 0.5 * rho * u * u;
    }
}

double CompressibleEulerSolver1D::pressureOf(const State& s, double gamma) {
    return (gamma - 1.0) * (s.E - 0.5 * s.rhoU * s.rhoU / s.rho);
}

double CompressibleEulerSolver1D::soundSpeedOf(const State& s, double gamma) {
    return std::sqrt(gamma * pressureOf(s, gamma) / s.rho);
}

CompressibleEulerSolver1D::State CompressibleEulerSolver1D::flux(const State& s, double gamma) {
    const double u = s.rhoU / s.rho;
    const double p = pressureOf(s, gamma);
    return State{s.rhoU, s.rhoU * u + p, u * (s.E + p)};
}

CompressibleEulerSolver1D::State CompressibleEulerSolver1D::rusanovFlux(const State& left,
                                                                         const State& right, double gamma) {
    const double uL = left.rhoU / left.rho;
    const double uR = right.rhoU / right.rho;
    const double sMax = std::max(std::fabs(uL) + soundSpeedOf(left, gamma),
                                  std::fabs(uR) + soundSpeedOf(right, gamma));
    const State fL = flux(left, gamma);
    const State fR = flux(right, gamma);
    return State{
        0.5 * (fL.rho + fR.rho) - 0.5 * sMax * (right.rho - left.rho),
        0.5 * (fL.rhoU + fR.rhoU) - 0.5 * sMax * (right.rhoU - left.rhoU),
        0.5 * (fL.E + fR.E) - 0.5 * sMax * (right.E - left.E),
    };
}

CompressibleEulerSolver1D::State CompressibleEulerSolver1D::stateAt(long long i) const {
    if (i >= 0 && i < static_cast<long long>(n_)) {
        const std::size_t idx = static_cast<std::size_t>(i);
        return State{rho_[idx], rhoU_[idx], E_[idx]};
    }
    // Outside the real domain: mirror the nearest real cell. Transmissive
    // copies it unchanged (an open end); Reflecting negates the momentum,
    // which forces u=0 exactly at that face by symmetry of the Rusanov
    // flux -- energy is untouched because kinetic energy (0.5*rho*u^2) does
    // not depend on the sign of u.
    const std::size_t neighbour = (i < 0) ? 0 : n_ - 1;
    const State interior{rho_[neighbour], rhoU_[neighbour], E_[neighbour]};
    if (boundary_ == BoundaryCondition::Transmissive) {
        return interior;
    }
    return State{interior.rho, -interior.rhoU, interior.E};
}

void CompressibleEulerSolver1D::step(double dt) {
    std::vector<double> newRho(n_);
    std::vector<double> newRhoU(n_);
    std::vector<double> newE(n_);

    // faceFlux starts as the left boundary face's flux and is carried
    // forward one face at a time, so each interior face is only evaluated
    // once rather than twice (as owner's right face and neighbour's left
    // face both).
    State faceFlux = rusanovFlux(stateAt(-1), stateAt(0), gamma_);
    for (std::size_t i = 0; i < n_; ++i) {
        const long long ii = static_cast<long long>(i);
        const State nextFlux = rusanovFlux(stateAt(ii), stateAt(ii + 1), gamma_);
        newRho[i] = rho_[i] - (dt / dx_) * (nextFlux.rho - faceFlux.rho);
        newRhoU[i] = rhoU_[i] - (dt / dx_) * (nextFlux.rhoU - faceFlux.rhoU);
        newE[i] = E_[i] - (dt / dx_) * (nextFlux.E - faceFlux.E);
        faceFlux = nextFlux;
    }

    rho_ = std::move(newRho);
    rhoU_ = std::move(newRhoU);
    E_ = std::move(newE);
    time_ += dt;
}

double CompressibleEulerSolver1D::stableTimeStep(double cfl) const {
    double maxSignalSpeed = 0.0;
    for (std::size_t i = 0; i < n_; ++i) {
        const State s{rho_[i], rhoU_[i], E_[i]};
        maxSignalSpeed = std::max(maxSignalSpeed, std::fabs(s.rhoU / s.rho) + soundSpeedOf(s, gamma_));
    }
    return cfl * dx_ / maxSignalSpeed;
}

double CompressibleEulerSolver1D::pressure(std::size_t cell) const {
    return pressureOf(State{rho_.at(cell), rhoU_.at(cell), E_.at(cell)}, gamma_);
}

double CompressibleEulerSolver1D::soundSpeed(std::size_t cell) const {
    return soundSpeedOf(State{rho_.at(cell), rhoU_.at(cell), E_.at(cell)}, gamma_);
}

double CompressibleEulerSolver1D::cellCenter(std::size_t cell) const {
    return (static_cast<double>(cell) + 0.5) * dx_;
}

double CompressibleEulerSolver1D::totalMass() const {
    double total = 0.0;
    for (double rho : rho_) {
        total += rho * dx_;
    }
    return total;
}

double CompressibleEulerSolver1D::totalMomentum() const {
    double total = 0.0;
    for (double rhoU : rhoU_) {
        total += rhoU * dx_;
    }
    return total;
}

double CompressibleEulerSolver1D::totalEnergy() const {
    double total = 0.0;
    for (double e : E_) {
        total += e * dx_;
    }
    return total;
}

} // namespace aether::solver
