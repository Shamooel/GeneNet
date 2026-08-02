#ifndef VARIANCECALCULATOR_H
#define VARIANCECALCULATOR_H

#include <vector>
#include "Gene.h"
#include "Sample.h"

// Class for manually calculating statistical metrics (mean & variance) for genes
class VarianceCalculator {
public:
    // Calculates mean expression for every gene across all samples
    static void calculateMeans(const std::vector<Sample>& samples, std::vector<Gene>& genes);

    // Calculates variance for every gene across all samples
    static void calculateVariances(const std::vector<Sample>& samples, std::vector<Gene>& genes);

    // Computes both mean and variance for all genes
    static void calculateAll(const std::vector<Sample>& samples, std::vector<Gene>& genes);
};

#endif // VARIANCECALCULATOR_H
