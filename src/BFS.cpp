#include "BFS.h"
#include <iostream>

using namespace std;

BFS::BFS() {}

void BFS::run(const Graph& graph) {
    componentSizes.clear();
    int numNodes = graph.getNumNodes();
    if (numNodes == 0) return;

    vector<bool> visited(numNodes, false);
    const vector<vector<Edge> >& adjList = graph.getAdjList();

    for (int i = 0; i < numNodes; ++i) {
        if (!visited[i]) {
            int nodeCount = 0;
            Queue q;

            visited[i] = true;
            q.enqueue(i);

            while (!q.isEmpty()) {
                int curr = q.front();
                q.dequeue();
                nodeCount++;

                const vector<Edge>& neighbors = adjList[curr];
                for (size_t k = 0; k < neighbors.size(); ++k) {
                    int neighborNode = neighbors[k].targetNode;
                    if (!visited[neighborNode]) {
                        visited[neighborNode] = true;
                        q.enqueue(neighborNode);
                    }
                }
            }
            componentSizes.push_back(nodeCount);
        }
    }
}

int BFS::getTotalComponents() const {
    return static_cast<int>(componentSizes.size());
}

const vector<int>& BFS::getComponentSizes() const {
    return componentSizes;
}

void BFS::displayResults() const {
    cout << "--------------------------------\n";
    cout << "BFS Completed\n";
    cout << "--------------------------------\n";
    cout << "Connected Components :\n\n";

    for (size_t i = 0; i < componentSizes.size(); ++i) {
        cout << "Component " << (i + 1) << " : " << componentSizes[i] << " Nodes\n";
    }
    cout << "\n";
}
