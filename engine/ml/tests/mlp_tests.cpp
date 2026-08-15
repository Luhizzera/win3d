#include "aether/mesh/StructuredGrid3D.hpp"
#include "aether/ml/MultiLayerPerceptron.hpp"
#include "aether/solver/SteadyDiffusionSolver.hpp"
#include "aether/testing/Check.hpp"

#include <cmath>
#include <cstdio>

using aether::core::Vector3;
using aether::mesh::StructuredGrid3D;
using aether::ml::MultiLayerPerceptron;
using aether::solver::SteadyDiffusionSolver;

namespace {

// The gold-standard correctness proof for a hand-rolled backprop
// implementation: compare its analytic gradient against a central finite
// difference of the exact same loss() function it is supposed to be the
// derivative of. Purely internal -- no external reference, no recalled
// fact, just calculus applied to a function this file can already
// evaluate directly. A bug in the chain-rule derivation (wrong sign, wrong
// activation derivative, wrong layer indexing) would show up here as a
// large mismatch on essentially every parameter, not a subtle one.
void testGradientMatchesFiniteDifference() {
    MultiLayerPerceptron net({2, 4, 1}, /*seed=*/7);
    const std::vector<std::vector<double>> inputs = {{0.3, -0.2}, {0.7, 0.1}, {-0.5, 0.4}, {0.0, 0.9}};
    const std::vector<std::vector<double>> targets = {{0.1}, {-0.3}, {0.6}, {0.2}};

    const std::vector<double> analytic = net.gradient(inputs, targets);
    const std::vector<double> original = net.parameters();

    const double eps = 1e-6;
    double maxDiff = 0.0;
    for (std::size_t i = 0; i < original.size(); ++i) {
        std::vector<double> plus = original;
        plus[i] += eps;
        net.setParameters(plus);
        const double lossPlus = net.loss(inputs, targets);

        std::vector<double> minus = original;
        minus[i] -= eps;
        net.setParameters(minus);
        const double lossMinus = net.loss(inputs, targets);

        net.setParameters(original);

        const double numeric = (lossPlus - lossMinus) / (2.0 * eps);
        maxDiff = std::max(maxDiff, std::fabs(numeric - analytic[i]));
    }
    std::printf("  [aether_ml_tests] max |analytic - finite-difference| gradient error = %.3e\n", maxDiff);
    AETHER_CHECK(maxDiff < 1e-4);
}

// XOR is the textbook proof that backprop + a hidden layer actually work,
// not just that the output layer's bias happens to fit the data: the four
// points (0,0)->0, (0,1)->1, (1,0)->1, (1,1)->0 are not linearly
// separable (no single line divides the two classes -- immediate from
// plotting the four points), so a network with no working hidden layer
// cannot solve this regardless of how its output weights are set. A bug
// that silently disabled the hidden layer's nonlinearity (e.g. an
// accidentally-identity activation) would pass a purely linear regression
// test but fail this one.
void testLearnsXor() {
    MultiLayerPerceptron net({2, 6, 1}, /*seed=*/3);
    const std::vector<std::vector<double>> inputs = {{0.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {1.0, 1.0}};
    const std::vector<std::vector<double>> targets = {{0.0}, {1.0}, {1.0}, {0.0}};

    // lr=0.5 was tried first and diverged to NaN by epoch 2000 -- plain
    // full-batch gradient descent overshooting with too large a step, not
    // a backprop bug (the gradient-check test above already proves the
    // gradient itself is correct to ~1e-11). Measured before picking a
    // final value, same as everywhere else in this project: lr=0.1 drives
    // the loss to exactly 0.000000 (5 decimal places) by epoch 2000 and
    // stays there, so 3000 epochs leaves a comfortable margin.
    double finalLoss = 0.0;
    for (int epoch = 0; epoch < 3000; ++epoch) {
        finalLoss = net.trainEpoch(inputs, targets, 0.1);
    }
    std::printf("  [aether_ml_tests] XOR training loss after 3000 epochs = %.6f\n", finalLoss);

    for (std::size_t i = 0; i < inputs.size(); ++i) {
        const double prediction = net.predict(inputs[i])[0];
        std::printf("  [aether_ml_tests] XOR(%.0f,%.0f) -> predicted=%.4f expected=%.0f\n", inputs[i][0],
                    inputs[i][1], prediction, targets[i][0]);
        AETHER_CHECK(std::fabs(prediction - targets[i][0]) < 0.05);
    }
}

// The actual point of Module 12's surrogate-model piece: train the network
// on real solver output, not a synthetic curve, and check it generalizes
// to inputs it never saw during training. Reuses the same
// SteadyDiffusionSolver "Poiseuille via source term" setup validated
// exactly in solver_tests.cpp and used by NelderMead's own inverse-problem
// test: value(15,0,0) = source * 1.05 (h=0.1, nx=30, target cell 15) --
// exactly linear in source, so a correctly-trained surrogate should
// reproduce it closely even on held-out points, and the check compares
// against both the real solver's own output *and* the closed-form
// prediction at those held-out points, not just one or the other.
double solvePoiseuilleValueAt15(double source) {
    const std::size_t nx = 30;
    const double h = 0.1;
    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(h * static_cast<double>(nx), 0.05, 0.05), nx, 1, 1);
    SteadyDiffusionSolver solver(grid);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, 0.0);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, 0.0);
    solver.setSourceTerm(source);
    solver.solveConjugateGradient(2000, 1e-12);
    return solver.value(15, 0, 0);
}

void testSurrogateGeneralizesToHeldOutSolverOutputs() {
    const double sourceMin = 0.5;
    const double sourceMax = 3.0;
    const double sourceMid = (sourceMin + sourceMax) / 2.0;
    const double sourceHalfRange = (sourceMax - sourceMin) / 2.0;
    const auto normalize = [&](double source) { return (source - sourceMid) / sourceHalfRange; };

    // A coarse training grid; the held-out test points below deliberately
    // fall *between* these training sources.
    std::vector<std::vector<double>> trainInputs;
    std::vector<std::vector<double>> trainTargets;
    for (double source = sourceMin; source <= sourceMax + 1e-9; source += 0.25) {
        const double value = solvePoiseuilleValueAt15(source);
        trainInputs.push_back({normalize(source)});
        trainTargets.push_back({value});
    }

    MultiLayerPerceptron net({1, 8, 1}, /*seed=*/11);
    double finalLoss = 0.0;
    for (int epoch = 0; epoch < 3000; ++epoch) {
        finalLoss = net.trainEpoch(trainInputs, trainTargets, 0.05);
    }
    std::printf("  [aether_ml_tests] Poiseuille surrogate training loss after 3000 epochs = %.6e\n", finalLoss);

    const std::vector<double> heldOutSources = {0.65, 1.4, 2.1, 2.85};
    for (double source : heldOutSources) {
        const double predicted = net.predict({normalize(source)})[0];
        const double actual = solvePoiseuilleValueAt15(source);
        const double closedForm = source * 1.05;
        AETHER_CHECK(std::fabs(actual - closedForm) < 1e-6); // sanity on the ground truth itself
        std::printf("  [aether_ml_tests] source=%.3f: surrogate=%.4f solver=%.4f closed-form=%.4f\n", source,
                    predicted, actual, closedForm);
        AETHER_CHECK(std::fabs(predicted - actual) < 0.05);
    }
}

} // namespace

int main() {
    testGradientMatchesFiniteDifference();
    testLearnsXor();
    testSurrogateGeneralizesToHeldOutSolverOutputs();
    std::printf("aether_ml_tests: OK\n");
    return 0;
}
