#ifndef PEARSONCORRELATION_H
#define PEARSONCORRELATION_H

#include <vector>
#include "Gene.h"
#include "Sample.h"

// Class for computing pairwise Pearson Correlation matrix for selected genes
class PearsonCorrelation {
private:
    std::vector<std::vector<double> > correlationMatrix;
    int totalPairsCompared;

public:
    PearsonCorrelation();

    // Computes single Pearson correlation coefficient between two gene index columns
    static double calculatePair(const std::vector<Sample>& samples,
                                int g1ColIdx, int g2ColIdx,
                                double mean1, double mean2,
                                double var1, double var2);

    // Computes correlation matrix for all selected genes
    void computeMatrix(const std::vector<Sample>& samples, const std::vector<Gene>& selectedGenes);

    // Returns correlation matrix reference
    const std::vector<std::vector<double> >& getMatrix() const;

    // Returns correlation value between selected gene i and selected gene j
    double getCorrelation(int i, int j) const;

    // Returns total unique gene pairs compared
    int getTotalPairsCompared() const;
};

#endif // PEARSONCORRELATION_H
