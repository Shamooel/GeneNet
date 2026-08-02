#ifndef BAGGING_H
#define BAGGING_H

#include <vector>
#include <string>
#include "Sample.h"
#include "Gene.h"
#include "DecisionTree.h"

// Class for Bagging ensemble classifier using bootstrap sampling, majority voting, and OOB Error calculation
class Bagging {
private:
    std::vector<DecisionTree*> trees;
    std::vector<std::vector<int> > bootstrapIndices;
    int numTrees;
    int treeMaxDepth;

    // Helper to generate a bootstrap sample with replacement
    std::vector<Sample> generateBootstrapSample(const std::vector<Sample>& samples,
                                                std::vector<int>& outIndices) const;

    // Helper to find majority vote among tree predictions
    std::string getMajorityVote(const std::vector<std::string>& votes) const;

public:
    Bagging(int nTrees = 5, int maxDepth = 5);
    ~Bagging();

    // Clears trained decision trees
    void clear();

    // Trains bagging ensemble
    void train(const std::vector<Sample>& samples, const std::vector<Gene>& selectedGenes);

    // Predicts sample class label using majority voting
    std::string predict(const Sample& sample) const;

    // Predicts class labels for a list of samples
    std::vector<std::string> predictAll(const std::vector<Sample>& samples) const;

    // Calculates Out-Of-Bag (OOB) Error percentage
    double calculateOOBError(const std::vector<Sample>& samples) const;

    // Traces detailed voting breakdown for a representative sample
    void traceBaggingPrediction(const Sample& sample,
                                const std::vector<std::string>& uniqueClasses) const;

    // Getters
    int getNumTrees() const;
    int getTreeMaxDepth() const;
};

#endif // BAGGING_H
