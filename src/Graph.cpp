#include "Graph.h"
#include <iostream>
#include <cmath>

using namespace std;

const double Graph::CORRELATION_THRESHOLD = 0.80;

Graph::Graph() : numNodes(0), numEdges(0) {}

void Graph::buildGraph(const vector<Gene>& selectedGenes,
                        const PearsonCorrelation& pearson,
                        double threshold) {
    nodes = selectedGenes;
    numNodes = static_cast<int>(nodes.size());
    numEdges = 0;

    adjList.clear();
    adjList.resize(numNodes);

    for (int i = 0; i < numNodes; ++i) {
        for (int j = i + 1; j < numNodes; ++j) {
            double r = pearson.getCorrelation(i, j);
            if (fabs(r) >= threshold) {
                adjList[i].push_back(Edge(j, r));
                adjList[j].push_back(Edge(i, r));
                numEdges++;
            }
        }
    }
}

int Graph::getNumNodes() const {
    return numNodes;
}

int Graph::getNumEdges() const {
    return numEdges;
}

double Graph::getAverageDegree() const {
    if (numNodes == 0) return 0.0;
    return (2.0 * static_cast<double>(numEdges)) / static_cast<double>(numNodes);
}

int Graph::getLargestDegree() const {
    int maxDeg = 0;
    for (size_t i = 0; i < adjList.size(); ++i) {
        int deg = static_cast<int>(adjList[i].size());
        if (deg > maxDeg) {
            maxDeg = deg;
        }
    }
    return maxDeg;
}

const vector<vector<Edge> >& Graph::getAdjList() const {
    return adjList;
}

const vector<Gene>& Graph::getNodes() const {
    return nodes;
}

void Graph::printGraphStats() const {
    cout << "--------------------------------\n";
    cout << "Graph Created Successfully\n";
    cout << "--------------------------------\n";
    cout << "Nodes: " << numNodes << "\n";
    cout << "Edges: " << numEdges << "\n";
    cout << "Average Degree: " << getAverageDegree() << "\n";
    cout << "Largest Degree: " << getLargestDegree() << "\n\n";
}
