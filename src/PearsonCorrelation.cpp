#include "PearsonCorrelation.h"
#include <cmath>

using namespace std;

PearsonCorrelation::PearsonCorrelation() : totalPairsCompared(0) {}

double PearsonCorrelation::calculatePair(const vector<Sample>& samples,
                                           int g1ColIdx, int g2ColIdx,
                                           double mean1, double mean2,
                                           double var1, double var2) {
    if (samples.empty()) return 0.0;

    double stdDevProduct = sqrt(var1) * sqrt(var2);
    if (stdDevProduct == 0.0) return 0.0;

    double covariance = 0.0;
    size_t numSamples = samples.size();

    for (size_t s = 0; s < numSamples; ++s) {
        const vector<double>& exprs = samples[s].getExpressions();
        double diff1 = exprs[g1ColIdx] - mean1;
        double diff2 = exprs[g2ColIdx] - mean2;
        covariance += diff1 * diff2;
    }

    covariance /= static_cast<double>(numSamples);
    return covariance / stdDevProduct;
}

void PearsonCorrelation::computeMatrix(const vector<Sample>& samples, const vector<Gene>& selectedGenes) {
    size_t K = selectedGenes.size();
    correlationMatrix.assign(K, vector<double>(K, 0.0));
    totalPairsCompared = 0;

    for (size_t i = 0; i < K; ++i) {
        correlationMatrix[i][i] = 1.0;
        for (size_t j = i + 1; j < K; ++j) {
            double r = calculatePair(samples,
                                      selectedGenes[i].getIndex(),
                                      selectedGenes[j].getIndex(),
                                      selectedGenes[i].getMean(),
                                      selectedGenes[j].getMean(),
                                      selectedGenes[i].getVariance(),
                                      selectedGenes[j].getVariance());
            correlationMatrix[i][j] = r;
            correlationMatrix[j][i] = r;
            totalPairsCompared++;
        }
    }
}

const vector<vector<double> >& PearsonCorrelation::getMatrix() const {
    return correlationMatrix;
}

double PearsonCorrelation::getCorrelation(int i, int j) const {
    if (i >= 0 && i < static_cast<int>(correlationMatrix.size()) &&
        j >= 0 && j < static_cast<int>(correlationMatrix.size())) {
        return correlationMatrix[i][j];
    }
    return 0.0;
}

int PearsonCorrelation::getTotalPairsCompared() const {
    return totalPairsCompared;
}
