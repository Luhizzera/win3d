#pragma once

#include <cstddef>
#include <vector>

namespace aether::ml {

// Module 12's "surrogate model" sub-area (first pass): a small feedforward
// multilayer perceptron for scalar/vector regression, trained by ordinary
// backpropagation + gradient descent, implemented from the algorithm's own
// definition rather than a wrapper over an external ML library. Same
// "hand-build the classic method from first principles" tradition as
// CG/GMRES/BiCGSTAB/Nelder-Mead elsewhere in this project -- unlike CUDA
// in Module 10, there is no real disadvantage to the dependency-free path
// here (a hand-rolled MLP costs some code, not throughput), so this
// doesn't reopen that kind of tradeoff.
//
// Intended use as a *surrogate*: train it to approximate a CFD solver's
// output as a function of some cheap input (e.g. a boundary condition or
// material parameter), so that afterwards `predict()` gives an answer
// without re-running the solver. The engine has no opinion on what the
// inputs/outputs mean -- same "caller decides the meaning" design as
// FieldArchive and NelderMead's objective function.
//
// Hidden layers use tanh; the output layer is linear (identity), the
// standard choice for regression rather than classification.
class MultiLayerPerceptron {
public:
    // layerSizes.front() is the input dimension, layerSizes.back() the
    // output dimension, everything between is a hidden layer's width.
    // Requires at least 2 entries (input and output, i.e. at least a
    // single-layer/no-hidden-layer network). Weights are initialized from
    // a seeded PRNG (uniform in +-1/sqrt(fanIn), the standard variance-
    // preserving scale) so construction is fully reproducible; biases
    // start at 0.
    explicit MultiLayerPerceptron(const std::vector<std::size_t>& layerSizes, unsigned seed = 42);

    std::vector<double> predict(const std::vector<double>& input) const;

    // Mean squared error (averaged over samples and output components)
    // over the given training set, at the network's current parameters.
    // Pure forward evaluation -- used both to report training progress and
    // as the reference value gradient() is checked against (finite
    // differences of this same function).
    double loss(const std::vector<std::vector<double>>& inputs, const std::vector<std::vector<double>>& targets) const;

    // Analytic gradient of loss() with respect to every weight and bias,
    // via backpropagation, flattened in the same order as parameters().
    std::vector<double> gradient(const std::vector<std::vector<double>>& inputs,
                                  const std::vector<std::vector<double>>& targets) const;

    // All weights and biases, flattened layer by layer (weights then
    // biases within each layer). Exists so gradient() and finite-
    // difference checks operate on the exact same representation.
    std::vector<double> parameters() const;
    void setParameters(const std::vector<double>& params);

    // One step of full-batch gradient descent. Returns loss() *before*
    // this step's update, so callers can watch it decrease across epochs.
    double trainEpoch(const std::vector<std::vector<double>>& inputs, const std::vector<std::vector<double>>& targets,
                       double learningRate);

private:
    struct ForwardCache {
        std::vector<std::vector<double>> activations;    // activations[0] = input, ..., activations[L] = output
        std::vector<std::vector<double>> preActivations;  // one per layer (pre-tanh/pre-identity)
    };
    ForwardCache forward(const std::vector<double>& input) const;

    std::vector<std::size_t> layerSizes_;
    std::vector<std::vector<double>> weights_; // per layer, row-major: weights_[l][o * inputs + i]
    std::vector<std::vector<double>> biases_;  // per layer
};

} // namespace aether::ml
