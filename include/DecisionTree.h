#ifndef DECISIONTREE_H
#define DECISIONTREE_H

#include <string>
#include <vector>
#include "Sample.h"
#include "Gene.h"

// Node structure for Decision Tree
struct TreeNode {
    int featureIndex;         // Gene index in sample's expression vector
    double threshold;         // Splitting threshold value
    std::string prediction;   // Class label prediction (for leaf nodes)
    TreeNode* left;           // Left child pointer (<= threshold)
    TreeNode* right;          // Right child pointer (> threshold)
    bool isLeaf;

    TreeNode()
        : featureIndex(-1), threshold(0.0), prediction(""),
          left(NULL), right(NULL), isLeaf(false) {}
};

// Class implementing a multi-class Decision Tree built with Gini Impurity
class DecisionTree {
private:
    TreeNode* root;
    int maxDepth;
    int minSamplesSplit;
    int internalNodes;
    int leafNodes;

    // Helper to free memory allocated for tree nodes
    void clear(TreeNode* node);

    // Calculates Gini Impurity for a set of sample pointers
    double calculateGini(const std::vector<const Sample*>& samplePtrs) const;

    // Finds most frequent class label among sample pointers
    std::string getMajorityClass(const std::vector<const Sample*>& samplePtrs) const;

    // Recursively constructs decision tree
    TreeNode* buildTree(const std::vector<const Sample*>& samplePtrs,
                        const std::vector<Gene>& selectedGenes,
                        int depth);

    // Helper for sample prediction traversal
    std::string predictSample(TreeNode* node, const Sample& sample) const;

    // Helper to print ASCII tree structure recursively
    void printNodeStructure(TreeNode* node,
                            const std::vector<Gene>& selectedGenes,
                            const std::string& prefix,
                            bool isLeft,
                            int currentDepth,
                            int maxDepthToPrint) const;

    // Helper to trace prediction path recursively
    void traceNodePath(TreeNode* node,
                       const Sample& sample,
                       const std::vector<Gene>& selectedGenes) const;

public:
    DecisionTree(int depth = 5, int minSamples = 2);
    ~DecisionTree();

    // Fits decision tree on samples using selected feature genes
    void train(const std::vector<Sample>& samples, const std::vector<Gene>& selectedGenes);

    // Predicts class label for a sample
    std::string predict(const Sample& sample) const;

    // Predicts class labels for multiple samples
    std::vector<std::string> predictAll(const std::vector<Sample>& samples) const;

    // Metric accessors
    int getMaxDepth() const;
    int getInternalNodes() const;
    int getLeafNodes() const;

    // Print first few levels of tree structure with ASCII branches
    void printTreeStructure(const std::vector<Gene>& selectedGenes, int maxDepthToPrint = 2) const;

    // Trace step-by-step decision path for a sample
    void tracePrediction(const Sample& sample, const std::vector<Gene>& selectedGenes) const;
};

#endif // DECISIONTREE_H
