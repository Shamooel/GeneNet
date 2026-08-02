#include "ProjectManager.h"
#include <iostream>
#include <sstream>

using namespace std;

ProjectManager::ProjectManager() {}

void ProjectManager::start() {
    int choice = 0;
    while (true) {
        menu.displayMainMenu();
        choice = menu.getChoice();

        if (choice == 1) {
            runPipeline();
        } else if (choice == 2) {
            runIndividualModule();
        } else if (choice == 3) {
            deleteResults();
        } else if (choice == 4) {
            cout << "\nExiting program...\n";
            break;
        } else {
            cout << "\nInvalid choice. Please try again.\n\n";
        }
    }
}

void ProjectManager::displayPipelineResults(const CSVReader& reader,
                                            const FeatureSelector& fs,
                                            const HashTable& ht,
                                            const Graph& g,
                                            const BFS& bfs,
                                            const DFS& dfs,
                                            const DecisionTree& dt,
                                            const Bagging& bagging,
                                            const Evaluation& eval) const {
    cout << "\n========================================\n";
    cout << "FULL PIPELINE RESULTS SUMMARY\n";
    cout << "========================================\n\n";

    cout << "--- Dataset Summary ---\n";
    cout << "Total Samples: " << reader.getSampleCount() << "\n";
    cout << "Total Genes: " << reader.getGeneCount() << "\n";
    cout << "Total Classes: " << reader.getClassCount() << "\n\n";

    fs.displayTopGenes(10);

    cout << "--- Hash Table Summary ---\n";
    cout << "Total Hash Table Entries: " << ht.getSize() << "\n\n";

    g.printGraphStats();

    dfs.displayResults();

    bfs.displayResults();

    cout << "--- Decision Tree ---\n";
    dt.printTreeStructure(fs.getSelectedGenes(), 2);
    if (!reader.getSamples().empty()) {
        // Pick a random representative sample so the trace varies per run
        int reprIdx = static_cast<int>(rand() % reader.getSamples().size());
        dt.tracePrediction(reader.getSamples()[reprIdx], fs.getSelectedGenes());

        cout << "--- Bagging Ensemble ---\n";
        bagging.traceBaggingPrediction(reader.getSamples()[reprIdx], reader.getUniqueClasses());
    } else {
        cout << "--- Bagging Ensemble ---\n";
    }

    cout << "--------------------------------\n";
    cout << "Evaluation Completed\n";
    cout << "--------------------------------\n";
    cout << "Note: Accuracy, Precision, and Recall are Training Set metrics.\n";
    cout << "Out-Of-Bag (OOB) Error provides the estimate of unseen-data performance.\n\n";
    cout << "Training Set Metrics:\n";
    cout << "  Accuracy:  " << eval.getAccuracy() << "%\n";
    cout << "  Precision: " << eval.getPrecision() << "%\n";
    cout << "  Recall:    " << eval.getRecall() << "%\n\n";
    cout << "Generalization Estimate (Unseen Data):\n";
    cout << "  Out-Of-Bag (OOB) Error Rate: " << bagging.calculateOOBError(reader.getSamples()) << "%\n\n";
}

void ProjectManager::runPipeline() {
    cout << "\nStep 1\nLoad Dataset\nRunning...\n";
    CSVReader csvReader;
    if (!csvReader.loadDataset("data/data.csv", "data/labels.csv")) {
        cout << "Error loading dataset.\n\n";
        return;
    }
    cout << "Completed\n\n";

    vector<Gene> genes = csvReader.getGenes();

    cout << "Step 2\nCalculate Mean\nRunning...\n";
    VarianceCalculator::calculateMeans(csvReader.getSamples(), genes);
    cout << "Completed\n\n";

    cout << "Step 3\nCalculate Variance\nRunning...\n";
    VarianceCalculator::calculateVariances(csvReader.getSamples(), genes);
    cout << "Completed\n\n";

    cout << "Step 4\nSort Genes\nRunning...\n";
    MergeSort::sort(genes);
    cout << "Completed\n\n";

    cout << "Step 5\nSelect Top Genes\nRunning...\n";
    FeatureSelector featureSelector;
    featureSelector.selectTopGenes(genes, 500);
    cout << "Completed\n\n";

    cout << "Step 6\nCreate Hash Table\nRunning...\n";
    HashTable geneHashTable;
    const vector<Gene>& selected = featureSelector.getSelectedGenes();
    for (size_t i = 0; i < selected.size(); ++i) {
        geneHashTable.insert(selected[i].getId(), static_cast<int>(i));
    }
    cout << "Completed\n\n";

    cout << "Step 7\nCompute Pearson Correlation\nRunning...\n";
    PearsonCorrelation pearson;
    pearson.computeMatrix(csvReader.getSamples(), selected);
    cout << "Completed\n\n";

    cout << "Step 8\nBuild Graph\nRunning...\n";
    Graph geneGraph;
    geneGraph.buildGraph(selected, pearson, Graph::CORRELATION_THRESHOLD);
    cout << "Completed\n\n";

    cout << "Step 9\nRun DFS\nRunning...\n";
    DFS dfs;
    dfs.run(geneGraph);
    cout << "Completed\n\n";

    cout << "Step 10\nRun BFS\nRunning...\n";
    BFS bfs;
    bfs.run(geneGraph);
    cout << "Completed\n\n";

    cout << "Step 11\nTrain Decision Tree\nRunning...\n";
    DecisionTree dt(5, 2);
    dt.train(csvReader.getSamples(), selected);
    cout << "Completed\n\n";

    cout << "Step 12\nTrain Bagging Ensemble\nRunning...\n";
    Bagging bagging(5, 5);
    bagging.train(csvReader.getSamples(), selected);
    cout << "Completed\n\n";

    cout << "Step 13\nEvaluate Performance\nRunning...\n";
    // Generate predictions using trained Bagging ensemble majority voting
    vector<string> predictions = bagging.predictAll(csvReader.getSamples());
    Evaluation eval;
    eval.evaluate(csvReader.getSamples(), predictions, csvReader.getUniqueClasses());
    cout << "Completed\n\n";

    cout << "====================================\n";
    cout << "Decision Tree\n";
    cout << "====================================\n";
    cout << "Training Completed\n\n";
    cout << "Maximum Depth : " << dt.getMaxDepth() << "\n";
    cout << "Internal Nodes : " << dt.getInternalNodes() << "\n";
    cout << "Leaf Nodes : " << dt.getLeafNodes() << "\n";
    cout << "Training Samples : " << csvReader.getSampleCount() << "\n";
    cout << "Selected Features : " << selected.size() << "\n";
    cout << "------------------------------------\n\n";

    cout << "====================================\n";
    cout << "Bagging\n";
    cout << "====================================\n";
    cout << "Bootstrap Samples Created\n";
    cout << "Decision Trees Trained\n";
    cout << "Majority Voting Completed\n";
    cout << "------------------------------------\n\n";

    cout << "--------------------------------\n";
    cout << "Evaluation Completed\n";
    cout << "--------------------------------\n";
    cout << "Note: Accuracy, Precision, and Recall are Training Set metrics.\n";
    cout << "Out-Of-Bag (OOB) Error provides the estimate of unseen-data performance.\n\n";
    cout << "Training Set Metrics:\n";
    cout << "  Accuracy:  " << eval.getAccuracy() << "%\n";
    cout << "  Precision: " << eval.getPrecision() << "%\n";
    cout << "  Recall:    " << eval.getRecall() << "%\n\n";
    cout << "Generalization Estimate (Unseen Data):\n";
    cout << "  Out-Of-Bag (OOB) Error Rate: " << bagging.calculateOOBError(csvReader.getSamples()) << "%\n\n";

    // Build text summary content to persist
    stringstream ss;
    ss << "Pipeline Run Execution Summary\n";
    ss << "Total Samples: " << csvReader.getSampleCount() << "\n";
    ss << "Total Genes: " << csvReader.getGeneCount() << "\n";
    ss << "Total Classes: " << csvReader.getClassCount() << "\n";
    ss << "Top Genes Selected: 500\n";
    ss << "Hash Table Entries: " << geneHashTable.getSize() << "\n";
    ss << "Pearson Pairs: " << pearson.getTotalPairsCompared() << "\n";
    ss << "Graph Nodes: " << geneGraph.getNumNodes() << "\n";
    ss << "Graph Edges: " << geneGraph.getNumEdges() << "\n";
    ss << "DFS Components: " << dfs.getTotalComponents() << "\n";
    ss << "BFS Components: " << bfs.getTotalComponents() << "\n";
    ss << "Training Accuracy: " << eval.getAccuracy() << "%\n";
    ss << "Training Precision: " << eval.getPrecision() << "%\n";
    ss << "Training Recall: " << eval.getRecall() << "%\n";
    ss << "OOB Error (Unseen Data Estimate): " << bagging.calculateOOBError(csvReader.getSamples()) << "%\n";

    string savedPath = fileManager.savePipelineRun(ss.str());
    cout << "Results saved to: " << savedPath << "\n";

    if (menu.askViewResults()) {
        displayPipelineResults(csvReader, featureSelector, geneHashTable, geneGraph, bfs, dfs, dt, bagging, eval);
    }
}

void ProjectManager::runIndividualModule() {
    while (true) {
        menu.displayIndividualModuleMenu();
        int modChoice = menu.getChoice();

        if (modChoice == 0) break;

        cout << "\nExecuting dependency chain for module selection...\n\n";

        // Step 1: Always load fresh dataset
        cout << "Step 1: Load Dataset\nRunning...\n";
        CSVReader csvReader;
        if (!csvReader.loadDataset("data/data.csv", "data/labels.csv")) {
            cout << "Error loading dataset.\n\n";
            return;
        }
        cout << "Completed\n\n";

        if (modChoice == 1) {
            stringstream ss;
            ss << "Module: CSV Reader\nSamples: " << csvReader.getSampleCount()
               << "\nGenes: " << csvReader.getGeneCount()
               << "\nClasses: " << csvReader.getClassCount() << "\n";
            string savedPath = fileManager.saveModuleRun("csv_reader", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                cout << "\nDataset Summary\n";
                cout << "Total Samples: " << csvReader.getSampleCount() << "\n";
                cout << "Total Genes: " << csvReader.getGeneCount() << "\n";
                cout << "Total Classes: " << csvReader.getClassCount() << "\n\n";
            }
            continue;
        }

        // Steps 2 & 3: Variance
        vector<Gene> genes = csvReader.getGenes();
        cout << "Step 2: Calculate Mean\nRunning...\n";
        VarianceCalculator::calculateMeans(csvReader.getSamples(), genes);
        cout << "Completed\n\n";

        cout << "Step 3: Calculate Variance\nRunning...\n";
        VarianceCalculator::calculateVariances(csvReader.getSamples(), genes);
        cout << "Completed\n\n";

        if (modChoice == 2) {
            stringstream ss;
            ss << "Module: Variance Calculation\nTotal Genes Processed: " << genes.size() << "\n";
            string savedPath = fileManager.saveModuleRun("variance", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                cout << "\nVariance Calculation Completed for " << genes.size() << " genes.\n\n";
            }
            continue;
        }

        // Steps 4 & 5: Sort & Feature Selection
        cout << "Step 4: Sort Genes\nRunning...\n";
        MergeSort::sort(genes);
        cout << "Completed\n\n";

        cout << "Step 5: Select Top Genes\nRunning...\n";
        FeatureSelector featureSelector;
        featureSelector.selectTopGenes(genes, 500);
        cout << "Completed\n\n";

        if (modChoice == 3) {
            stringstream ss;
            ss << "Module: Top Gene Selection\nTop 500 Genes Selected.\n";
            string savedPath = fileManager.saveModuleRun("feature_selection", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                featureSelector.displayTopGenes(10);
            }
            continue;
        }

        const vector<Gene>& selected = featureSelector.getSelectedGenes();

        if (modChoice == 4) { // Hash Table
            cout << "Step 6: Create Hash Table\nRunning...\n";
            HashTable ht;
            for (size_t i = 0; i < selected.size(); ++i) {
                ht.insert(selected[i].getId(), static_cast<int>(i));
            }
            cout << "Completed\n\n";

            stringstream ss;
            ss << "Module: Hash Table\nTotal Entries: " << ht.getSize() << "\n";
            string savedPath = fileManager.saveModuleRun("hash_table", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                cout << "\nHash Table Created\nTotal Entries: " << ht.getSize() << "\n\n";
            }
            continue;
        }

        if (modChoice == 5) { // Pearson Correlation
            cout << "Step 7: Compute Pearson Correlation\nRunning...\n";
            PearsonCorrelation pearson;
            pearson.computeMatrix(csvReader.getSamples(), selected);
            cout << "Completed\n\n";

            stringstream ss;
            ss << "Module: Pearson Correlation\nTotal Pairs Compared: " << pearson.getTotalPairsCompared() << "\n";
            string savedPath = fileManager.saveModuleRun("pearson", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                cout << "\nPearson Correlation Completed\nTotal Gene Pairs Compared: " << pearson.getTotalPairsCompared() << "\n\n";
            }
            continue;
        }

        // Dependencies for Graph, BFS, DFS
        cout << "Step 7: Compute Pearson Correlation\nRunning...\n";
        PearsonCorrelation pearson;
        pearson.computeMatrix(csvReader.getSamples(), selected);
        cout << "Completed\n\n";

        cout << "Step 8: Build Graph\nRunning...\n";
        Graph geneGraph;
        geneGraph.buildGraph(selected, pearson, Graph::CORRELATION_THRESHOLD);
        cout << "Completed\n\n";

        if (modChoice == 6) { // Graph
            stringstream ss;
            ss << "Module: Graph\nNodes: " << geneGraph.getNumNodes()
               << "\nEdges: " << geneGraph.getNumEdges()
               << "\nAvg Degree: " << geneGraph.getAverageDegree()
               << "\nLargest Degree: " << geneGraph.getLargestDegree() << "\n";
            string savedPath = fileManager.saveModuleRun("graph", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                geneGraph.printGraphStats();
            }
            continue;
        }

        if (modChoice == 7) { // BFS
            cout << "Step 9: Run BFS\nRunning...\n";
            BFS bfs;
            bfs.run(geneGraph);
            cout << "Completed\n\n";

            stringstream ss;
            ss << "Module: BFS\nConnected Components: " << bfs.getTotalComponents() << "\n";
            string savedPath = fileManager.saveModuleRun("bfs", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                bfs.displayResults();
            }
            continue;
        }

        if (modChoice == 8) { // DFS
            cout << "Step 9: Run DFS\nRunning...\n";
            DFS dfs;
            dfs.run(geneGraph);
            cout << "Completed\n\n";

            stringstream ss;
            ss << "Module: DFS\nConnected Components: " << dfs.getTotalComponents() << "\n";
            string savedPath = fileManager.saveModuleRun("dfs", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                dfs.displayResults();
            }
            continue;
        }

        if (modChoice == 9) { // Decision Tree
            cout << "Step 11: Decision Tree\nRunning...\n";
            DecisionTree dt(5, 2);
            dt.train(csvReader.getSamples(), selected);
            cout << "Completed\n\n";

            stringstream ss;
            ss << "Module: Decision Tree\nTrained Depth 5 Decision Tree.\n";
            string savedPath = fileManager.saveModuleRun("decision_tree", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                dt.printTreeStructure(selected, 2);
                if (!csvReader.getSamples().empty()) {
                    int reprIdx = static_cast<int>(rand() % csvReader.getSamples().size());
                    dt.tracePrediction(csvReader.getSamples()[reprIdx], selected);
                }
            }
            continue;
        }

        if (modChoice == 10) { // Bagging
            cout << "Step 12: Bagging\nRunning...\n";
            Bagging bagging(5, 5);
            bagging.train(csvReader.getSamples(), selected);
            cout << "Completed\n\n";

            stringstream ss;
            ss << "Module: Bagging\nTrained 5 Trees Bagging Ensemble.\n";
            string savedPath = fileManager.saveModuleRun("bagging", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                if (!csvReader.getSamples().empty()) {
                    int reprIdx = static_cast<int>(rand() % csvReader.getSamples().size());
                    bagging.traceBaggingPrediction(csvReader.getSamples()[reprIdx], csvReader.getUniqueClasses());
                }
            }
            continue;
        }

        if (modChoice == 11) { // Evaluation
            cout << "Step 12: Bagging\nRunning...\n";
            Bagging bagging(5, 5);
            bagging.train(csvReader.getSamples(), selected);
            cout << "Completed\n\n";

            cout << "Step 13: Evaluation\nRunning...\n";
            vector<string> predictions = bagging.predictAll(csvReader.getSamples());
            Evaluation eval;
            eval.evaluate(csvReader.getSamples(), predictions, csvReader.getUniqueClasses());
            cout << "Completed\n\n";

            stringstream ss;
            ss << "Module: Evaluation\nAccuracy (Training): " << eval.getAccuracy() << "%\n"
               << "Precision (Training): " << eval.getPrecision() << "%\n"
               << "Recall (Training): " << eval.getRecall() << "%\n";
            string savedPath = fileManager.saveModuleRun("evaluation", ss.str());
            cout << "Module result saved to: " << savedPath << "\n";

            if (menu.askViewResults()) {
                eval.displayResults();
            }
            continue;
        }

        cout << "Invalid module choice.\n\n";
    }
}

void ProjectManager::deleteResults() {
    while (true) {
        menu.displayDeleteMenu();
        int delChoice = menu.getChoice();

        if (delChoice == 0) break;

        if (delChoice == 1) {
            fileManager.deletePipelineResults();
        } else if (delChoice == 2) {
            fileManager.deleteModuleResults();
        } else if (delChoice == 3) {
            fileManager.deleteAllResults();
        } else {
            cout << "Invalid choice. Please try again.\n\n";
        }
    }
}
