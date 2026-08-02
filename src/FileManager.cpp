#include "FileManager.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

FileManager::FileManager() {}

void FileManager::createDirectory(const string& path) const {
    string cmd = "mkdir \"" + path + "\" 2>NUL";
    system(cmd.c_str());
}

void FileManager::removeDirectory(const string& path) const {
    string cmd = "rmdir /S /Q \"" + path + "\" 2>NUL";
    system(cmd.c_str());
}

bool FileManager::confirmDeletion(const string& message) const {
    cout << message << " (y/n): ";
    char choice;
    cin >> choice;
    if (choice == 'y' || choice == 'Y') {
        return true;
    }
    return false;
}

int FileManager::getNextPipelineRunNumber() const {
    int run = 1;
    while (true) {
        stringstream ss;
        ss << "output/pipeline/pipeline_run_" << run << "/results.txt";
        ifstream test(ss.str().c_str());
        if (!test.good()) {
            break;
        }
        test.close();
        run++;
    }
    return run;
}

int FileManager::getNextModuleRunNumber(const string& moduleName) const {
    int run = 1;
    while (true) {
        stringstream ss;
        ss << "output/modules/" << moduleName << "/run_" << run << ".txt";
        ifstream test(ss.str().c_str());
        if (!test.good()) {
            break;
        }
        test.close();
        run++;
    }
    return run;
}

string FileManager::savePipelineRun(const string& content) {
    createDirectory("output");
    createDirectory("output/pipeline");

    int runNum = getNextPipelineRunNumber();
    stringstream dirSS;
    dirSS << "output/pipeline/pipeline_run_" << runNum;
    string runDir = dirSS.str();
    createDirectory(runDir);

    string filePath = runDir + "/results.txt";
    ofstream outFile(filePath.c_str());
    if (outFile.is_open()) {
        outFile << content;
        outFile.close();
    }
    return filePath;
}

string FileManager::saveModuleRun(const string& moduleName, const string& content) {
    createDirectory("output");
    createDirectory("output/modules");

    string modDir = "output/modules/" + moduleName;
    createDirectory(modDir);

    int runNum = getNextModuleRunNumber(moduleName);
    stringstream ss;
    ss << modDir << "/run_" << runNum << ".txt";
    string filePath = ss.str();

    ofstream outFile(filePath.c_str());
    if (outFile.is_open()) {
        outFile << content;
        outFile.close();
    }
    return filePath;
}

void FileManager::deletePipelineResults() {
    if (confirmDeletion("Are you sure you want to delete all Pipeline results?")) {
        removeDirectory("output\\pipeline");
        cout << "Pipeline results deleted successfully.\n\n";
    } else {
        cout << "Deletion canceled.\n\n";
    }
}

void FileManager::deleteModuleResults() {
    if (confirmDeletion("Are you sure you want to delete all Individual Module results?")) {
        removeDirectory("output\\modules");
        cout << "Module results deleted successfully.\n\n";
    } else {
        cout << "Deletion canceled.\n\n";
    }
}

void FileManager::deleteAllResults() {
    if (confirmDeletion("Are you sure you want to delete ALL output results?")) {
        removeDirectory("output");
        cout << "All results deleted successfully.\n\n";
    } else {
        cout << "Deletion canceled.\n\n";
    }
}
