#include "Bagging.h"
#include <iostream>
#include <cstdlib>

using namespace std;

Bagging::Bagging(int nTrees, int maxDepth)
    : numTrees(nTrees), treeMaxDepth(maxDepth) {}

Bagging::~Bagging() {
    clear();
}

void Bagging::clear() {
    for (size_t i = 0; i < trees.size(); ++i) {
        delete trees[i];
    }
    trees.clear();
    bootstrapIndices.clear();
}

vector<Sample> Bagging::generateBootstrapSample(const vector<Sample>& samples,
                                                vector<int>& outIndices) const {
    vector<Sample> bootstrap;
    size_t n = samples.size();
    bootstrap.reserve(n);
    outIndices.clear();
    outIndices.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        int randIdx = rand() % n;
        bootstrap.push_back(samples[randIdx]);
        outIndices.push_back(randIdx);
    }
    return bootstrap;
}

string Bagging::getMajorityVote(const vector<string>& votes) const {
    if (votes.empty()) return "";

    vector<string> labels;
    vector<int> counts;

    for (size_t i = 0; i < votes.size(); ++i) {
        string v = votes[i];
        int foundIdx = -1;
        for (size_t j = 0; j < labels.size(); ++j) {
            if (labels[j] == v) {
                foundIdx = static_cast<int>(j);
                break;
            }
        }
        if (foundIdx != -1) {
            counts[foundIdx]++;
        } else {
            labels.push_back(v);
            counts.push_back(1);
        }
    }

    int maxCount = -1;
    string winner = "";
    for (size_t i = 0; i < labels.size(); ++i) {
        if (counts[i] > maxCount) {
            maxCount = counts[i];
            winner = labels[i];
        }
    }
    return winner;
}

void Bagging::train(const vector<Sample>& samples, const vector<Gene>& selectedGenes) {
    clear();
    if (samples.empty()) return;

    // RNG is seeded once at program startup in main(); do NOT reseed here.
    for (int t = 0; t < numTrees; ++t) {
        vector<int> sampleIdxs;
        vector<Sample> bootstrapSample = generateBootstrapSample(samples, sampleIdxs);
        bootstrapIndices.push_back(sampleIdxs);

        DecisionTree* tree = new DecisionTree(treeMaxDepth, 2);
        tree->train(bootstrapSample, selectedGenes);
        trees.push_back(tree);
    }
}

string Bagging::predict(const Sample& sample) const {
    if (trees.empty()) return "";

    vector<string> votes;
    votes.reserve(trees.size());

    for (size_t i = 0; i < trees.size(); ++i) {
        votes.push_back(trees[i]->predict(sample));
    }
    return getMajorityVote(votes);
}

vector<string> Bagging::predictAll(const vector<Sample>& samples) const {
    vector<string> predictions;
    predictions.reserve(samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
        predictions.push_back(predict(samples[i]));
    }
    return predictions;
}

double Bagging::calculateOOBError(const vector<Sample>& samples) const {
    if (trees.empty() || samples.empty()) return 0.0;

    int oobTested = 0;
    int oobCorrect = 0;

    for (size_t i = 0; i < samples.size(); ++i) {
        vector<string> oobVotes;

        for (size_t t = 0; t < trees.size(); ++t) {
            bool inBootstrap = false;
            for (size_t b = 0; b < bootstrapIndices[t].size(); ++b) {
                if (bootstrapIndices[t][b] == static_cast<int>(i)) {
                    inBootstrap = true;
                    break;
                }
            }
            if (!inBootstrap) {
                oobVotes.push_back(trees[t]->predict(samples[i]));
            }
        }

        if (!oobVotes.empty()) {
            oobTested++;
            string oobPred = getMajorityVote(oobVotes);
            if (oobPred == samples[i].getLabel()) {
                oobCorrect++;
            }
        }
    }

    if (oobTested == 0) return 0.0;
    double oobAcc = (static_cast<double>(oobCorrect) / static_cast<double>(oobTested)) * 100.0;
    return 100.0 - oobAcc;
}

void Bagging::traceBaggingPrediction(const Sample& sample,
                                     const vector<string>& uniqueClasses) const {
    cout << "====================================\n";
    cout << "Bagging Prediction Example\n";
    cout << "====================================\n";
    cout << "Sample ID : " << sample.getId() << "\n";
    cout << "Actual Class : " << sample.getLabel() << "\n\n";

    vector<string> votes;
    for (size_t i = 0; i < trees.size(); ++i) {
        string pred = trees[i]->predict(sample);
        votes.push_back(pred);
        cout << "Tree " << (i + 1) << " -> " << pred << "\n";
    }

    cout << "\n------------------------------------\n";
    cout << "Votes\n";
    for (size_t c = 0; c < uniqueClasses.size(); ++c) {
        string cls = uniqueClasses[c];
        int count = 0;
        for (size_t v = 0; v < votes.size(); ++v) {
            if (votes[v] == cls) {
                count++;
            }
        }
        cout << cls << " : " << count << "\n";
    }

    string finalPred = getMajorityVote(votes);
    cout << "------------------------------------\n";
    cout << "Final Prediction : " << finalPred << "\n";
    cout << "Result : " << (finalPred == sample.getLabel() ? "Correct" : "Incorrect") << "\n\n";
}

int Bagging::getNumTrees() const {
    return numTrees;
}

int Bagging::getTreeMaxDepth() const {
    return treeMaxDepth;
}
