#include "DecisionTree.h"
#include <iostream>

using namespace std;

DecisionTree::DecisionTree(int depth, int minSamples)
    : root(NULL), maxDepth(depth), minSamplesSplit(minSamples),
      internalNodes(0), leafNodes(0) {}

DecisionTree::~DecisionTree() {
    clear(root);
}

void DecisionTree::clear(TreeNode* node) {
    if (node == NULL) return;
    clear(node->left);
    clear(node->right);
    delete node;
}

static int getClassIndex(const string& label) {
    if (label == "PRAD") return 0;
    if (label == "LUAD") return 1;
    if (label == "BRCA") return 2;
    if (label == "KIRC") return 3;
    if (label == "COAD") return 4;
    return 5;
}

string DecisionTree::getMajorityClass(const vector<const Sample*>& samplePtrs) const {
    if (samplePtrs.empty()) return "";

    vector<string> labels;
    vector<int> counts;

    for (size_t i = 0; i < samplePtrs.size(); ++i) {
        string lbl = samplePtrs[i]->getLabel();
        int foundIdx = -1;
        for (size_t j = 0; j < labels.size(); ++j) {
            if (labels[j] == lbl) {
                foundIdx = static_cast<int>(j);
                break;
            }
        }
        if (foundIdx != -1) {
            counts[foundIdx]++;
        } else {
            labels.push_back(lbl);
            counts.push_back(1);
        }
    }

    int maxCount = -1;
    string bestLabel = "";
    for (size_t i = 0; i < labels.size(); ++i) {
        if (counts[i] > maxCount) {
            maxCount = counts[i];
            bestLabel = labels[i];
        }
    }
    return bestLabel;
}

double DecisionTree::calculateGini(const vector<const Sample*>& samplePtrs) const {
    if (samplePtrs.empty()) return 0.0;

    int counts[8] = {0};
    for (size_t i = 0; i < samplePtrs.size(); ++i) {
        int idx = getClassIndex(samplePtrs[i]->getLabel());
        counts[idx]++;
    }

    double total = static_cast<double>(samplePtrs.size());
    double sumSqProb = 0.0;
    for (int c = 0; c < 8; ++c) {
        if (counts[c] > 0) {
            double p = static_cast<double>(counts[c]) / total;
            sumSqProb += (p * p);
        }
    }
    return 1.0 - sumSqProb;
}

TreeNode* DecisionTree::buildTree(const vector<const Sample*>& samplePtrs,
                                  const vector<Gene>& selectedGenes,
                                  int depth) {
    TreeNode* node = new TreeNode();

    double currentGini = calculateGini(samplePtrs);
    string majorityLabel = getMajorityClass(samplePtrs);
    node->prediction = majorityLabel;

    // Stopping conditions: pure node, max depth, or min samples
    if (currentGini == 0.0 || depth >= maxDepth || static_cast<int>(samplePtrs.size()) < minSamplesSplit) {
        node->isLeaf = true;
        leafNodes++;
        return node;
    }

    double bestGiniGain = -1.0;
    int bestGeneColIdx = -1;
    double bestThreshold = 0.0;
    vector<const Sample*> bestLeftPtrs;
    vector<const Sample*> bestRightPtrs;

    size_t numFeatures = selectedGenes.size();
    size_t numSamplesInNode = samplePtrs.size();

    // For each feature, compute 5 statistically-spread threshold candidates via a
    // single linear scan (min, 25%, 50%, 75%, max-eps).  This is O(n) per feature,
    // uses no prohibited containers, and ensures different bootstrap samples produce
    // different candidates because bootstrap shifts the min/max/median of each feature.
    static const double FRACTIONS[5] = { 0.0, 0.25, 0.5, 0.75, 1.0 };

    for (size_t f = 0; f < numFeatures; ++f) {
        int colIdx = selectedGenes[f].getIndex();

        // Single linear scan to find min and max for this feature in the current node
        double minVal = samplePtrs[0]->getExpressions()[colIdx];
        double maxVal = minVal;
        for (size_t i = 1; i < numSamplesInNode; ++i) {
            double v = samplePtrs[i]->getExpressions()[colIdx];
            if (v < minVal) minVal = v;
            if (v > maxVal) maxVal = v;
        }

        // Skip constant features — all values identical, no useful split exists
        if (maxVal - minVal < 1e-12) continue;

        double range = maxVal - minVal;

        for (int ci = 0; ci < 5; ++ci) {
            // Candidate threshold: fraction of [min, max-eps] so the right child
            // is always non-empty when threshold == maxVal
            double threshold = minVal + FRACTIONS[ci] * range;
            if (ci == 4) threshold = maxVal - 1e-9; // just below max

            vector<const Sample*> leftPtrs;
            vector<const Sample*> rightPtrs;

            for (size_t i = 0; i < numSamplesInNode; ++i) {
                if (samplePtrs[i]->getExpressions()[colIdx] <= threshold) {
                    leftPtrs.push_back(samplePtrs[i]);
                } else {
                    rightPtrs.push_back(samplePtrs[i]);
                }
            }

            if (leftPtrs.empty() || rightPtrs.empty()) continue;

            double leftGini  = calculateGini(leftPtrs);
            double rightGini = calculateGini(rightPtrs);
            double splitGini = (static_cast<double>(leftPtrs.size())  / numSamplesInNode) * leftGini +
                               (static_cast<double>(rightPtrs.size()) / numSamplesInNode) * rightGini;

            double giniGain = currentGini - splitGini;

            // Tie-breaking: when Gini Gain is equal (within floating-point epsilon),
            // prefer the more balanced split.  Balance is measured by the size of the
            // smaller child — a larger smaller-child means a more even partition.
            // Different bootstrap samples produce different balances, so ties break
            // differently across trees, naturally increasing ensemble diversity.
            size_t thisBalance  = leftPtrs.size() < rightPtrs.size() ? leftPtrs.size() : rightPtrs.size();
            size_t bestBalance  = bestLeftPtrs.size() < bestRightPtrs.size() ? bestLeftPtrs.size() : bestRightPtrs.size();

            bool strictlyBetter = giniGain > bestGiniGain + 1e-12;
            bool tieBreakBetter = (giniGain >= bestGiniGain - 1e-12) &&
                                  (giniGain <= bestGiniGain + 1e-12) &&
                                  (thisBalance > bestBalance);

            if (strictlyBetter || tieBreakBetter) {
                bestGiniGain   = giniGain;
                bestGeneColIdx = colIdx;
                bestThreshold  = threshold;
                bestLeftPtrs   = leftPtrs;
                bestRightPtrs  = rightPtrs;
            }
        }
    }

    if (bestGiniGain <= 0.0 || bestLeftPtrs.empty() || bestRightPtrs.empty()) {
        node->isLeaf = true;
        leafNodes++;
        return node;
    }

    node->featureIndex = bestGeneColIdx;
    node->threshold = bestThreshold;
    node->isLeaf = false;
    internalNodes++;

    node->left = buildTree(bestLeftPtrs, selectedGenes, depth + 1);
    node->right = buildTree(bestRightPtrs, selectedGenes, depth + 1);

    return node;
}

void DecisionTree::train(const vector<Sample>& samples, const vector<Gene>& selectedGenes) {
    clear(root);
    root = NULL;
    internalNodes = 0;
    leafNodes = 0;

    vector<const Sample*> samplePtrs;
    samplePtrs.reserve(samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
        samplePtrs.push_back(&samples[i]);
    }
    root = buildTree(samplePtrs, selectedGenes, 0);
}

string DecisionTree::predictSample(TreeNode* node, const Sample& sample) const {
    if (node == NULL) return "";
    if (node->isLeaf || node->left == NULL || node->right == NULL) {
        return node->prediction;
    }

    if (sample.getExpressions()[node->featureIndex] <= node->threshold) {
        return predictSample(node->left, sample);
    } else {
        return predictSample(node->right, sample);
    }
}

string DecisionTree::predict(const Sample& sample) const {
    return predictSample(root, sample);
}

vector<string> DecisionTree::predictAll(const vector<Sample>& samples) const {
    vector<string> predictions;
    predictions.reserve(samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
        predictions.push_back(predict(samples[i]));
    }
    return predictions;
}

int DecisionTree::getMaxDepth() const {
    return maxDepth;
}

int DecisionTree::getInternalNodes() const {
    return internalNodes;
}

int DecisionTree::getLeafNodes() const {
    return leafNodes;
}

// Find gene name for a column index
static string findGeneName(int colIdx, const vector<Gene>& selectedGenes) {
    for (size_t i = 0; i < selectedGenes.size(); ++i) {
        if (selectedGenes[i].getIndex() == colIdx) {
            return selectedGenes[i].getId();
        }
    }
    return "Gene_" + to_string(colIdx);
}

void DecisionTree::printNodeStructure(TreeNode* node,
                                      const vector<Gene>& selectedGenes,
                                      const string& prefix,
                                      bool isLeft,
                                      int currentDepth,
                                      int maxDepthToPrint) const {
    if (node == NULL || currentDepth > maxDepthToPrint) return;

    if (currentDepth == 0) {
        cout << "Root\n";
        if (!node->isLeaf) {
            cout << "  Gene : " << findGeneName(node->featureIndex, selectedGenes) << "\n";
            cout << "  Threshold : " << node->threshold << "\n";
        } else {
            cout << "  Predict : " << node->prediction << "\n";
        }
    } else {
        // ASCII branch characters — compatible with Dev-C++ 5.11 console (CP850/1252)
        string branch = isLeft ? "+-- Left" : "`-- Right";
        cout << prefix << branch << "\n";
        if (node->isLeaf) {
            cout << prefix << "  Predict : " << node->prediction << "\n";
        } else {
            cout << prefix << "  Gene : " << findGeneName(node->featureIndex, selectedGenes) << "\n";
            cout << prefix << "  Threshold : " << node->threshold << "\n";
        }
    }

    if (!node->isLeaf && currentDepth < maxDepthToPrint) {
        // Use plain '|' instead of Unicode box-drawing character
        string nextPrefix = prefix + (isLeft ? "|   " : "    ");
        if (currentDepth == 0) nextPrefix = "";
        printNodeStructure(node->left, selectedGenes, nextPrefix, true, currentDepth + 1, maxDepthToPrint);
        printNodeStructure(node->right, selectedGenes, nextPrefix, false, currentDepth + 1, maxDepthToPrint);
    }
}

void DecisionTree::printTreeStructure(const vector<Gene>& selectedGenes, int maxDepthToPrint) const {
    cout << "====================================\n";
    cout << "Decision Tree Structure\n";
    cout << "====================================\n";
    if (root == NULL) {
        cout << "Tree not trained.\n\n";
        return;
    }
    printNodeStructure(root, selectedGenes, "", true, 0, maxDepthToPrint);
    cout << "\n";
}

void DecisionTree::traceNodePath(TreeNode* node,
                                 const Sample& sample,
                                 const vector<Gene>& selectedGenes) const {
    if (node == NULL) return;
    if (node->isLeaf || node->left == NULL || node->right == NULL) {
        cout << "\n";
        cout << "Predicted Class : " << node->prediction << "\n";
        cout << "Actual Class    : " << sample.getLabel() << "\n";
        cout << "Result          : " << (node->prediction == sample.getLabel() ? "Correct" : "Incorrect") << "\n\n";
        return;
    }

    double val = sample.getExpressions()[node->featureIndex];
    string geneName = findGeneName(node->featureIndex, selectedGenes);

    // ASCII trace: plain comparison lines + explicit direction label
    cout << "Gene " << geneName << " = " << val << "\n";
    if (val <= node->threshold) {
        cout << val << " <= " << node->threshold << "\n";
        cout << "Move to LEFT child\n";
        cout << "-->\n";
        traceNodePath(node->left, sample, selectedGenes);
    } else {
        cout << val << " > " << node->threshold << "\n";
        cout << "Move to RIGHT child\n";
        cout << "-->\n";
        traceNodePath(node->right, sample, selectedGenes);
    }
}

void DecisionTree::tracePrediction(const Sample& sample, const vector<Gene>& selectedGenes) const {
    cout << "====================================\n";
    cout << "Decision Tree Prediction Example\n";
    cout << "====================================\n";
    cout << "Sample ID : " << sample.getId() << "\n";
    traceNodePath(root, sample, selectedGenes);
}
