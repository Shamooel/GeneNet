#ifndef PROJECTMANAGER_H
#define PROJECTMANAGER_H

#include "Menu.h"
#include "CSVReader.h"
#include "VarianceCalculator.h"
#include "MergeSort.h"
#include "FeatureSelector.h"
#include "HashTable.h"
#include "PearsonCorrelation.h"
#include "Graph.h"
#include "BFS.h"
#include "DFS.h"
#include "DecisionTree.h"
#include "Bagging.h"
#include "Evaluation.h"
#include "FileManager.h"

// Class controlling program flow, menu handling, dependency chain execution, and outputs
class ProjectManager {
private:
    Menu menu;
    FileManager fileManager;

    // Execution methods
    void runPipeline();
    void runIndividualModule();
    void deleteResults();

    // Helper to display full pipeline result summary
    void displayPipelineResults(const CSVReader& reader,
                                const FeatureSelector& fs,
                                const HashTable& ht,
                                const Graph& g,
                                const BFS& bfs,
                                const DFS& dfs,
                                const DecisionTree& dt,
                                const Bagging& bagging,
                                const Evaluation& eval) const;

public:
    ProjectManager();

    // Main entry point for application execution loop
    void start();
};

#endif // PROJECTMANAGER_H
