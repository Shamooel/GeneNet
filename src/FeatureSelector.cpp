#include "FeatureSelector.h"
#include <iostream>

using namespace std;

FeatureSelector::FeatureSelector() {}

void FeatureSelector::selectTopGenes(const vector<Gene>& sortedGenes, int topK) {
    selectedGenes.clear();
    int limit = topK;
    if (limit > static_cast<int>(sortedGenes.size())) {
        limit = static_cast<int>(sortedGenes.size());
    }

    for (int i = 0; i < limit; ++i) {
        selectedGenes.push_back(sortedGenes[i]);
    }
}

const vector<Gene>& FeatureSelector::getSelectedGenes() const {
    return selectedGenes;
}

void FeatureSelector::displayTopGenes(int count) const {
    cout << "--------------------------------\n";
    cout << "Top " << count << " Most Variable Genes\n";
    cout << "--------------------------------\n";
    cout << "Rank\tGene Name\tVariance\tOriginal Column Index\n";

    int limit = count;
    if (limit > static_cast<int>(selectedGenes.size())) {
        limit = static_cast<int>(selectedGenes.size());
    }

    for (int i = 0; i < limit; ++i) {
        cout << (i + 1) << "\t"
             << selectedGenes[i].getId() << "\t\t"
             << selectedGenes[i].getVariance() << "\t\t"
             << selectedGenes[i].getIndex() << "\n";
    }
    cout << "\n";
}
