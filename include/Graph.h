#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include "Gene.h"
#include "PearsonCorrelation.h"
#include "HashTable.h"

// Struct representing an Edge in an Adjacency List
struct Edge {
    int targetNode;
    double weight; // Pearson correlation coefficient

    Edge() : targetNode(-1), weight(0.0) {}
    Edge(int target, double w) : targetNode(target), weight(w) {}
};

// Gene Co-expression Graph implemented purely with an Adjacency List
class Graph {
private:
    int numNodes;
    int numEdges;
    std::vector<std::vector<Edge> > adjList; // ADJACENCY LIST ONLY
    std::vector<Gene> nodes;

public:
    static const double CORRELATION_THRESHOLD; // 0.80

    Graph();

    // Builds dynamic Adjacency List graph from selected genes and Pearson correlations
    void buildGraph(const std::vector<Gene>& selectedGenes,
                    const PearsonCorrelation& pearson,
                    double threshold = 0.80);

    // Accessors for graph statistics
    int getNumNodes() const;
    int getNumEdges() const;
    double getAverageDegree() const;
    int getLargestDegree() const;

    // Returns reference to Adjacency List
    const std::vector<std::vector<Edge> >& getAdjList() const;

    // Returns list of node genes
    const std::vector<Gene>& getNodes() const;

    // Displays formatted graph output statistics
    void printGraphStats() const;
};

#endif // GRAPH_H
