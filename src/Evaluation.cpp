#include "Evaluation.h"
#include <iostream>

using namespace std;

Evaluation::Evaluation() : accuracy(0.0), precision(0.0), recall(0.0) {}

void Evaluation::evaluate(const vector<Sample>& samples,
                          const vector<string>& predictions,
                          const vector<string>& uniqueClasses) {
    if (samples.empty() || samples.size() != predictions.size()) return;

    size_t totalSamples = samples.size();
    int correctCount = 0;

    for (size_t i = 0; i < totalSamples; ++i) {
        if (samples[i].getLabel() == predictions[i]) {
            correctCount++;
        }
    }

    accuracy = (static_cast<double>(correctCount) / static_cast<double>(totalSamples)) * 100.0;

    double sumPrecision = 0.0;
    double sumRecall = 0.0;
    int validClassesCount = 0;

    for (size_t c = 0; c < uniqueClasses.size(); ++c) {
        string targetClass = uniqueClasses[c];
        int tp = 0;
        int fp = 0;
        int fn = 0;

        for (size_t i = 0; i < totalSamples; ++i) {
            string actual = samples[i].getLabel();
            string pred = predictions[i];

            if (actual == targetClass && pred == targetClass) {
                tp++;
            } else if (actual != targetClass && pred == targetClass) {
                fp++;
            } else if (actual == targetClass && pred != targetClass) {
                fn++;
            }
        }

        double classPrecision = 0.0;
        if (tp + fp > 0) {
            classPrecision = static_cast<double>(tp) / static_cast<double>(tp + fp);
        }

        double classRecall = 0.0;
        if (tp + fn > 0) {
            classRecall = static_cast<double>(tp) / static_cast<double>(tp + fn);
        }

        sumPrecision += classPrecision;
        sumRecall += classRecall;
        validClassesCount++;
    }

    if (validClassesCount > 0) {
        precision = (sumPrecision / static_cast<double>(validClassesCount)) * 100.0;
        recall = (sumRecall / static_cast<double>(validClassesCount)) * 100.0;
    } else {
        precision = 0.0;
        recall = 0.0;
    }
}

double Evaluation::getAccuracy() const {
    return accuracy;
}

double Evaluation::getPrecision() const {
    return precision;
}

double Evaluation::getRecall() const {
    return recall;
}

void Evaluation::displayResults() const {
    cout << "--------------------------------\n";
    cout << "Evaluation Completed\n";
    cout << "--------------------------------\n";
    cout << "Accuracy: " << accuracy << "%\n";
    cout << "Precision: " << precision << "%\n";
    cout << "Recall: " << recall << "%\n\n";
}
