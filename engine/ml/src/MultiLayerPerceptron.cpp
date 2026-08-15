#include "aether/ml/MultiLayerPerceptron.hpp"

#include <cmath>
#include <random>
#include <stdexcept>

namespace aether::ml {

namespace {

double tanhDerivativeFromOutput(double tanhValue) { return 1.0 - tanhValue * tanhValue; }

} // namespace

MultiLayerPerceptron::MultiLayerPerceptron(const std::vector<std::size_t>& layerSizes, unsigned seed)
    : layerSizes_(layerSizes) {
    if (layerSizes_.size() < 2) {
        throw std::invalid_argument("MultiLayerPerceptron: layerSizes must have at least an input and output layer");
    }

    std::mt19937 rng(seed);
    const std::size_t numLayers = layerSizes_.size() - 1;
    weights_.resize(numLayers);
    biases_.resize(numLayers);
    for (std::size_t l = 0; l < numLayers; ++l) {
        const std::size_t fanIn = layerSizes_[l];
        const std::size_t fanOut = layerSizes_[l + 1];
        const double limit = 1.0 / std::sqrt(static_cast<double>(fanIn));
        std::uniform_real_distribution<double> dist(-limit, limit);

        weights_[l].resize(fanOut * fanIn);
        for (double& w : weights_[l]) {
            w = dist(rng);
        }
        biases_[l].assign(fanOut, 0.0);
    }
}

MultiLayerPerceptron::ForwardCache MultiLayerPerceptron::forward(const std::vector<double>& input) const {
    if (input.size() != layerSizes_.front()) {
        throw std::invalid_argument("MultiLayerPerceptron: input size does not match the network's input layer");
    }

    const std::size_t numLayers = weights_.size();
    ForwardCache cache;
    cache.activations.resize(numLayers + 1);
    cache.preActivations.resize(numLayers);
    cache.activations[0] = input;

    for (std::size_t l = 0; l < numLayers; ++l) {
        const std::size_t fanIn = layerSizes_[l];
        const std::size_t fanOut = layerSizes_[l + 1];
        const std::vector<double>& prevActivation = cache.activations[l];

        std::vector<double> z(fanOut);
        for (std::size_t o = 0; o < fanOut; ++o) {
            double sum = biases_[l][o];
            for (std::size_t i = 0; i < fanIn; ++i) {
                sum += weights_[l][o * fanIn + i] * prevActivation[i];
            }
            z[o] = sum;
        }
        cache.preActivations[l] = z;

        std::vector<double> a(fanOut);
        const bool isOutputLayer = (l + 1 == numLayers);
        for (std::size_t o = 0; o < fanOut; ++o) {
            a[o] = isOutputLayer ? z[o] : std::tanh(z[o]);
        }
        cache.activations[l + 1] = a;
    }

    return cache;
}

std::vector<double> MultiLayerPerceptron::predict(const std::vector<double>& input) const {
    return forward(input).activations.back();
}

double MultiLayerPerceptron::loss(const std::vector<std::vector<double>>& inputs,
                                   const std::vector<std::vector<double>>& targets) const {
    if (inputs.size() != targets.size() || inputs.empty()) {
        throw std::invalid_argument("MultiLayerPerceptron::loss: inputs/targets must be non-empty and same size");
    }
    const std::size_t outputDim = layerSizes_.back();

    double sumSquaredError = 0.0;
    for (std::size_t n = 0; n < inputs.size(); ++n) {
        const std::vector<double> prediction = predict(inputs[n]);
        for (std::size_t o = 0; o < outputDim; ++o) {
            const double diff = prediction[o] - targets[n][o];
            sumSquaredError += diff * diff;
        }
    }
    return sumSquaredError / static_cast<double>(inputs.size() * outputDim);
}

std::vector<double> MultiLayerPerceptron::gradient(const std::vector<std::vector<double>>& inputs,
                                                     const std::vector<std::vector<double>>& targets) const {
    if (inputs.size() != targets.size() || inputs.empty()) {
        throw std::invalid_argument("MultiLayerPerceptron::gradient: inputs/targets must be non-empty and same size");
    }
    const std::size_t numLayers = weights_.size();
    const std::size_t outputDim = layerSizes_.back();
    const double normalization = 1.0 / static_cast<double>(inputs.size() * outputDim);

    std::vector<std::vector<double>> weightGrad(numLayers);
    std::vector<std::vector<double>> biasGrad(numLayers);
    for (std::size_t l = 0; l < numLayers; ++l) {
        weightGrad[l].assign(weights_[l].size(), 0.0);
        biasGrad[l].assign(biases_[l].size(), 0.0);
    }

    for (std::size_t n = 0; n < inputs.size(); ++n) {
        const ForwardCache cache = forward(inputs[n]);

        // delta = dL/d(pre-activation) for the current layer, starting at
        // the output layer (linear activation, so d(pre-activation) =
        // d(activation) there).
        std::vector<double> delta(outputDim);
        for (std::size_t o = 0; o < outputDim; ++o) {
            delta[o] = 2.0 * (cache.activations[numLayers][o] - targets[n][o]) * normalization;
        }

        for (std::size_t layerIndex = numLayers; layerIndex-- > 0;) {
            const std::size_t fanIn = layerSizes_[layerIndex];
            const std::size_t fanOut = layerSizes_[layerIndex + 1];
            const std::vector<double>& prevActivation = cache.activations[layerIndex];

            for (std::size_t o = 0; o < fanOut; ++o) {
                biasGrad[layerIndex][o] += delta[o];
                for (std::size_t i = 0; i < fanIn; ++i) {
                    weightGrad[layerIndex][o * fanIn + i] += delta[o] * prevActivation[i];
                }
            }

            if (layerIndex == 0) {
                break; // no earlier layer to propagate into
            }

            std::vector<double> deltaPrev(fanIn, 0.0);
            for (std::size_t i = 0; i < fanIn; ++i) {
                double sum = 0.0;
                for (std::size_t o = 0; o < fanOut; ++o) {
                    sum += weights_[layerIndex][o * fanIn + i] * delta[o];
                }
                // prevActivation feeds into a tanh layer (every layer
                // except the output is tanh), so multiply by that
                // activation's own derivative before it becomes the next
                // layer's delta.
                deltaPrev[i] = sum * tanhDerivativeFromOutput(prevActivation[i]);
            }
            delta = std::move(deltaPrev);
        }
    }

    std::vector<double> flat;
    for (std::size_t l = 0; l < numLayers; ++l) {
        flat.insert(flat.end(), weightGrad[l].begin(), weightGrad[l].end());
        flat.insert(flat.end(), biasGrad[l].begin(), biasGrad[l].end());
    }
    return flat;
}

std::vector<double> MultiLayerPerceptron::parameters() const {
    std::vector<double> flat;
    for (std::size_t l = 0; l < weights_.size(); ++l) {
        flat.insert(flat.end(), weights_[l].begin(), weights_[l].end());
        flat.insert(flat.end(), biases_[l].begin(), biases_[l].end());
    }
    return flat;
}

void MultiLayerPerceptron::setParameters(const std::vector<double>& params) {
    std::size_t offset = 0;
    for (std::size_t l = 0; l < weights_.size(); ++l) {
        for (double& w : weights_[l]) {
            w = params.at(offset++);
        }
        for (double& b : biases_[l]) {
            b = params.at(offset++);
        }
    }
    if (offset != params.size()) {
        throw std::invalid_argument("MultiLayerPerceptron::setParameters: parameter count does not match the network");
    }
}

double MultiLayerPerceptron::trainEpoch(const std::vector<std::vector<double>>& inputs,
                                         const std::vector<std::vector<double>>& targets, double learningRate) {
    const double lossBefore = loss(inputs, targets);
    const std::vector<double> grad = gradient(inputs, targets);
    std::vector<double> params = parameters();
    for (std::size_t i = 0; i < params.size(); ++i) {
        params[i] -= learningRate * grad[i];
    }
    setParameters(params);
    return lossBefore;
}

} // namespace aether::ml
