#ifndef FEATURESELECTOR_H
#define FEATURESELECTOR_H

#include <vector>
#include "Gene.h"

// Class for feature selection (extracting top variable genes)
class FeatureSelector {
private:
    std::vector<Gene> selectedGenes;

public:
    FeatureSelector();

    // Selects Top K (default 500) genes from sorted list
    void selectTopGenes(const std::vector<Gene>& sortedGenes, int topK = 500);

    // Returns vector of selected Top K genes
    const std::vector<Gene>& getSelectedGenes() const;

    // Displays formatted table of Top N (default 10) most variable genes
    void displayTopGenes(int count = 10) const;
};

#endif // FEATURESELECTOR_H
