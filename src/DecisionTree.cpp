#include "DecisionTree.h"
#include <iostream>
#include <cmath>

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

    // Dynamic class counting via stack arrays — no hardcoded labels, no prohibited
    // containers, supports any number of cancer classes up to 32.
    string classLabels[32];
    int    classCounts[32];
    int    numClasses = 0;

    for (size_t i = 0; i < samplePtrs.size(); ++i) {
        const string& lbl = samplePtrs[i]->getLabel();
        int found = -1;
        for (int c = 0; c < numClasses; ++c) {
            if (classLabels[c] == lbl) { found = c; break; }
        }
        if (found >= 0) {
            classCounts[found]++;
        } else if (numClasses < 32) {
            classLabels[numClasses] = lbl;
            classCounts[numClasses] = 1;
            numClasses++;
        }
    }

    double total = static_cast<double>(samplePtrs.size());
    double sumSqProb = 0.0;
    for (int c = 0; c < numClasses; ++c) {
        double p = static_cast<double>(classCounts[c]) / total;
        sumSqProb += p * p;
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
    // Use epsilon comparison: FP accumulation can yield ~1e-16 instead of exactly 0.0
    if (currentGini < 1e-12 || depth >= maxDepth || static_cast<int>(samplePtrs.size()) < minSamplesSplit) {
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

    // For each feature, compute 9 threshold candidates using a single O(n) linear
    // scan to collect min, max, sum, and sum-of-squares.  This yields:
    //   coarse grid  : 0%, 10%, 20%, 30%
    //   distribution : mean-sigma, mean, mean+sigma  (data-aware, varies per bootstrap)
    //   tails        : 80%, max-eps
    // No std::sort. No prohibited containers. All stack-allocated.

    for (size_t f = 0; f < numFeatures; ++f) {
        int colIdx = selectedGenes[f].getIndex();

        // Single pass: collect min, max, sum, sumSq
        double minVal = samplePtrs[0]->getExpressions()[colIdx];
        double maxVal = minVal;
        double sum    = minVal;
        double sumSq  = minVal * minVal;

        for (size_t i = 1; i < numSamplesInNode; ++i) {
            double v = samplePtrs[i]->getExpressions()[colIdx];
            if (v < minVal) minVal = v;
            if (v > maxVal) maxVal = v;
            sum   += v;
            sumSq += v * v;
        }

        // Skip constant features: every sample has the same value, no split is useful
        if (maxVal - minVal < 1e-12) continue;

        double range = maxVal - minVal;
        double mean  = sum / static_cast<double>(numSamplesInNode);
        double var   = (sumSq / static_cast<double>(numSamplesInNode)) - (mean * mean);
        if (var < 0.0) var = 0.0;          // guard against tiny FP-negative variance
        double sigma = sqrt(var);

        // Build 9 candidates on the stack
        double candidates[9];
        candidates[0] = minVal;
        candidates[1] = minVal + range * 0.1;
        candidates[2] = minVal + range * 0.2;
        candidates[3] = minVal + range * 0.3;
        candidates[4] = (sigma > 1e-12) ? (mean - sigma) : (minVal + range * 0.4);
        candidates[5] = mean;
        candidates[6] = (sigma > 1e-12) ? (mean + sigma) : (minVal + range * 0.6);
        candidates[7] = minVal + range * 0.8;
        candidates[8] = maxVal - 1e-9;     // just below max: right child always non-empty

        for (int ci = 0; ci < 9; ++ci) {
            double threshold = candidates[ci];
            // Clamp: sigma extremes can push candidates outside [min, max)
            if (threshold < minVal)           threshold = minVal;
            if (threshold >= maxVal - 1e-9)   threshold = maxVal - 1e-9;

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

            // Tie-breaking: prefer more balanced split when gains are within epsilon.
            // Bootstrap changes each tree's distribution, so balances differ per tree,
            // naturally breaking ties in different directions across the ensemble.
            size_t thisBalance = leftPtrs.size() < rightPtrs.size() ? leftPtrs.size() : rightPtrs.size();
            size_t bestBalance = bestLeftPtrs.size() < bestRightPtrs.size() ? bestLeftPtrs.size() : bestRightPtrs.size();

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

    if (bestGiniGain < 1e-12 || bestLeftPtrs.empty() || bestRightPtrs.empty()) {
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

// Returns the name of the gene used at the root split
string DecisionTree::getRootFeatureName(const vector<Gene>& selectedGenes) const {
    if (root == NULL || root->isLeaf) return "N/A";
    return findGeneName(root->featureIndex, selectedGenes);
}

// Returns the threshold value at the root split
double DecisionTree::getRootThreshold() const {
    if (root == NULL || root->isLeaf) return 0.0;
    return root->threshold;
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
        // Build child prefix: root's direct children start with ""; deeper nodes
        // inherit the parent prefix plus a branch-continuation character.
        string nextPrefix;
        if (currentDepth == 0) {
            nextPrefix = "";          // root's children: no leading indent
        } else {
            nextPrefix = prefix + (isLeft ? "|   " : "    ");
        }
        printNodeStructure(node->left,  selectedGenes, nextPrefix, true,  currentDepth + 1, maxDepthToPrint);
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
