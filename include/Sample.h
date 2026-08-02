#ifndef SAMPLE_H
#define SAMPLE_H

#include <string>
#include <vector>

// Class representing a biological Sample with ID, cancer label, and gene expressions
class Sample {
private:
    std::string id;
    std::string label;
    std::vector<double> expressions;

public:
    // Default constructor
    Sample() : id(""), label("") {}

    // Parameterized constructor
    Sample(const std::string& sampleId, const std::string& sampleLabel = "")
        : id(sampleId), label(sampleLabel) {}

    // Getter for Sample ID
    std::string getId() const {
        return id;
    }

    // Getter for Cancer Label
    std::string getLabel() const {
        return label;
    }

    // Getter for expression values
    const std::vector<double>& getExpressions() const {
        return expressions;
    }

    // Setter for Sample ID
    void setId(const std::string& sampleId) {
        id = sampleId;
    }

    // Setter for Cancer Label
    void setLabel(const std::string& sampleLabel) {
        label = sampleLabel;
    }

    // Setter for full expression vector
    void setExpressions(const std::vector<double>& exprs) {
        expressions = exprs;
    }

    // Add a single expression value
    void addExpression(double val) {
        expressions.push_back(val);
    }
};

#endif // SAMPLE_H
