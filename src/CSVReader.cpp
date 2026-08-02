#include "CSVReader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>

using namespace std;

CSVReader::CSVReader() {}

bool CSVReader::isClassPresent(const string& className) const {
    for (size_t i = 0; i < uniqueClasses.size(); ++i) {
        if (uniqueClasses[i] == className) {
            return true;
        }
    }
    return false;
}

bool CSVReader::parseLabels(const string& labelsFilePath) {
    ifstream file(labelsFilePath.c_str());
    if (!file.is_open()) {
        cerr << "Error: Could not open labels file: " << labelsFilePath << endl;
        return false;
    }

    samples.clear();
    uniqueClasses.clear();

    string line;
    if (!getline(file, line)) {
        file.close();
        return false;
    }

    while (getline(file, line)) {
        if (line.empty()) continue;

        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) continue;

        size_t commaPos = line.find(',');
        if (commaPos == string::npos) continue;

        string sampleId = line.substr(0, commaPos);
        string label = line.substr(commaPos + 1);

        Sample s(sampleId, label);
        samples.push_back(s);

        if (!label.empty() && !isClassPresent(label)) {
            uniqueClasses.push_back(label);
        }
    }

    file.close();
    return true;
}

bool CSVReader::parseData(const string& dataFilePath) {
    ifstream file(dataFilePath.c_str());
    if (!file.is_open()) {
        cerr << "Error: Could not open data file: " << dataFilePath << endl;
        return false;
    }

    genes.clear();

    string line;
    if (!getline(file, line)) {
        file.close();
        return false;
    }

    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }

    stringstream headerStream(line);
    string token;

    getline(headerStream, token, ',');

    int geneIdx = 0;
    while (getline(headerStream, token, ',')) {
        genes.push_back(Gene(token, geneIdx));
        geneIdx++;
    }

    size_t sampleRowIdx = 0;
    while (getline(file, line)) {
        if (line.empty()) continue;

        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) continue;

        const char* str = line.c_str();
        const char* comma = strchr(str, ',');
        if (!comma) continue;

        string sampleId(str, comma - str);

        vector<double> exprs;
        exprs.reserve(genes.size());

        const char* ptr = comma + 1;
        char* endPtr = nullptr;
        while (*ptr != '\0') {
            double val = strtod(ptr, &endPtr);
            if (endPtr == ptr) {
                if (*ptr == ',') {
                    ptr++;
                    continue;
                }
                break;
            }
            exprs.push_back(val);
            ptr = endPtr;
            if (*ptr == ',') {
                ptr++;
            }
        }

        if (sampleRowIdx < samples.size()) {
            samples[sampleRowIdx].setExpressions(exprs);
        } else {
            Sample s(sampleId, "");
            s.setExpressions(exprs);
            samples.push_back(s);
        }
        sampleRowIdx++;
    }

    file.close();
    return true;
}

bool CSVReader::loadDataset(const string& dataFilePath, const string& labelsFilePath) {
    if (!parseLabels(labelsFilePath)) {
        return false;
    }
    if (!parseData(dataFilePath)) {
        return false;
    }
    return true;
}

int CSVReader::getSampleCount() const {
    return static_cast<int>(samples.size());
}

int CSVReader::getGeneCount() const {
    return static_cast<int>(genes.size());
}

int CSVReader::getClassCount() const {
    return static_cast<int>(uniqueClasses.size());
}

const vector<Sample>& CSVReader::getSamples() const {
    return samples;
}

const vector<Gene>& CSVReader::getGenes() const {
    return genes;
}

const vector<string>& CSVReader::getUniqueClasses() const {
    return uniqueClasses;
}
