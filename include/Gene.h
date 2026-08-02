#ifndef GENE_H
#define GENE_H

#include <string>

// Class representing a single Gene entity with ID, column index, mean, and variance
class Gene {
private:
    std::string id;
    int index;
    double mean;
    double variance;

public:
    // Default constructor
    Gene() : id(""), index(-1), mean(0.0), variance(0.0) {}

    // Parameterized constructor
    Gene(const std::string& geneId, int geneIndex)
        : id(geneId), index(geneIndex), mean(0.0), variance(0.0) {}

    // Full constructor
    Gene(const std::string& geneId, int geneIndex, double geneMean, double geneVar)
        : id(geneId), index(geneIndex), mean(geneMean), variance(geneVar) {}

    // Getter for Gene ID
    std::string getId() const {
        return id;
    }

    // Getter for original column index
    int getIndex() const {
        return index;
    }

    // Getter for Mean
    double getMean() const {
        return mean;
    }

    // Getter for Variance
    double getVariance() const {
        return variance;
    }

    // Setter for Gene ID
    void setId(const std::string& geneId) {
        id = geneId;
    }

    // Setter for column index
    void setIndex(int geneIndex) {
        index = geneIndex;
    }

    // Setter for Mean
    void setMean(double geneMean) {
        mean = geneMean;
    }

    // Setter for Variance
    void setVariance(double geneVar) {
        variance = geneVar;
    }
};

#endif // GENE_H
