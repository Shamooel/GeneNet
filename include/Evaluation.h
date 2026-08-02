#ifndef EVALUATION_H
#define EVALUATION_H

#include <vector>
#include <string>
#include "Sample.h"

// Class for computing classification performance metrics (Accuracy, Precision, Recall)
class Evaluation {
private:
    double accuracy;
    double precision;
    double recall;

public:
    Evaluation();

    // Evaluates actual vs predicted labels across unique cancer classes
    void evaluate(const std::vector<Sample>& samples,
                  const std::vector<std::string>& predictions,
                  const std::vector<std::string>& uniqueClasses);

    // Get calculated Accuracy (in percentage: 0 to 100)
    double getAccuracy() const;

    // Get calculated Macro Precision (in percentage: 0 to 100)
    double getPrecision() const;

    // Get calculated Macro Recall (in percentage: 0 to 100)
    double getRecall() const;

    // Displays formatted evaluation metrics
    void displayResults() const;
};

#endif // EVALUATION_H
