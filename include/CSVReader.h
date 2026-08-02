#ifndef CSVREADER_H
#define CSVREADER_H

#include <string>
#include <vector>
#include "Gene.h"
#include "Sample.h"

// Class responsible for reading dataset CSV files and extracting samples, genes, and classes
class CSVReader {
private:
    std::vector<Sample> samples;
    std::vector<Gene> genes;
    std::vector<std::string> uniqueClasses;

    // Helper method to parse labels CSV file
    bool parseLabels(const std::string& labelsFilePath);

    // Helper method to parse expression data CSV file
    bool parseData(const std::string& dataFilePath);

    // Helper method to check if a class string is already in uniqueClasses list without using std::find
    bool isClassPresent(const std::string& className) const;

public:
    CSVReader();

    // Main method to load full dataset
    bool loadDataset(const std::string& dataFilePath, const std::string& labelsFilePath);

    // Accessors for calculated statistics
    int getSampleCount() const;
    int getGeneCount() const;
    int getClassCount() const;

    // Accessors for parsed data
    const std::vector<Sample>& getSamples() const;
    const std::vector<Gene>& getGenes() const;
    const std::vector<std::string>& getUniqueClasses() const;
};

#endif // CSVREADER_H
