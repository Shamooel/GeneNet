#include "VarianceCalculator.h"
#include <cmath>

using namespace std;

void VarianceCalculator::calculateMeans(const vector<Sample>& samples, vector<Gene>& genes) {
    if (samples.empty() || genes.empty()) return;

    size_t numSamples = samples.size();
    size_t numGenes = genes.size();

    for (size_t g = 0; g < numGenes; ++g) {
        double sum = 0.0;
        for (size_t s = 0; s < numSamples; ++s) {
            const vector<double>& exprs = samples[s].getExpressions();
            if (g < exprs.size()) {
                sum += exprs[g];
            }
        }
        double meanVal = sum / static_cast<double>(numSamples);
        genes[g].setMean(meanVal);
    }
}

void VarianceCalculator::calculateVariances(const vector<Sample>& samples, vector<Gene>& genes) {
    if (samples.empty() || genes.empty()) return;

    size_t numSamples = samples.size();
    size_t numGenes = genes.size();

    for (size_t g = 0; g < numGenes; ++g) {
        double geneMean = genes[g].getMean();
        double sumSquaredDiff = 0.0;

        for (size_t s = 0; s < numSamples; ++s) {
            const vector<double>& exprs = samples[s].getExpressions();
            if (g < exprs.size()) {
                double diff = exprs[g] - geneMean;
                sumSquaredDiff += diff * diff;
            }
        }
        double varianceVal = sumSquaredDiff / static_cast<double>(numSamples);
        genes[g].setVariance(varianceVal);
    }
}

void VarianceCalculator::calculateAll(const vector<Sample>& samples, vector<Gene>& genes) {
    calculateMeans(samples, genes);
    calculateVariances(samples, genes);
}
