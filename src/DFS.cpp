#include "DFS.h"
#include <iostream>

using namespace std;

DFS::DFS() {}

void DFS::run(const Graph& graph) {
    componentSizes.clear();
    int numNodes = graph.getNumNodes();
    if (numNodes == 0) return;

    vector<bool> visited(numNodes, false);
    const vector<vector<Edge> >& adjList = graph.getAdjList();

    for (int i = 0; i < numNodes; ++i) {
        if (!visited[i]) {
            int nodeCount = 0;
            Stack s;

            s.push(i);

            while (!s.isEmpty()) {
                int curr = s.top();
                s.pop();

                if (!visited[curr]) {
                    visited[curr] = true;
                    nodeCount++;

                    const vector<Edge>& neighbors = adjList[curr];
                    for (size_t k = 0; k < neighbors.size(); ++k) {
                        int neighborNode = neighbors[k].targetNode;
                        if (!visited[neighborNode]) {
                            s.push(neighborNode);
                        }
                    }
                }
            }
            componentSizes.push_back(nodeCount);
        }
    }
}

int DFS::getTotalComponents() const {
    return static_cast<int>(componentSizes.size());
}

const vector<int>& DFS::getComponentSizes() const {
    return componentSizes;
}

void DFS::displayResults() const {
    cout << "--------------------------------\n";
    cout << "DFS Completed\n";
    cout << "--------------------------------\n";
    cout << "Connected Components :\n\n";

    for (size_t i = 0; i < componentSizes.size(); ++i) {
        cout << "Component " << (i + 1) << " : " << componentSizes[i] << " Nodes\n";
    }
    cout << "\n";
}
