#ifndef FILEMANAGER_H
#define FILEMANAGER_H

#include <string>

// Class managing file system persistence for pipeline and module execution runs
class FileManager {
private:
    // Creates directory safely using system command
    void createDirectory(const std::string& path) const;

    // Determines next available run number for complete pipeline
    int getNextPipelineRunNumber() const;

    // Determines next available run number for an individual module
    int getNextModuleRunNumber(const std::string& moduleName) const;

    // Removes directory recursively
    void removeDirectory(const std::string& path) const;

    // Prompts user for (y/n) deletion confirmation
    bool confirmDeletion(const std::string& message) const;

public:
    FileManager();

    // Saves complete pipeline output summary to output/pipeline/pipeline_run_X/
    std::string savePipelineRun(const std::string& content);

    // Saves individual module output summary to output/modules/<moduleName>/run_X.txt
    std::string saveModuleRun(const std::string& moduleName, const std::string& content);

    // Deletes all complete pipeline results after user confirmation
    void deletePipelineResults();

    // Deletes all individual module results after user confirmation
    void deleteModuleResults();

    // Deletes all output results after user confirmation
    void deleteAllResults();
};

#endif // FILEMANAGER_H
