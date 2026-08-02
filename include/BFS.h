#ifndef BFS_H
#define BFS_H

#include <vector>
#include "Graph.h"
#include "Queue.h"

// Class for Breadth-First Search connected component traversal using custom Queue
class BFS {
private:
    std::vector<int> componentSizes;

public:
    BFS();

    // Runs BFS traversal on graph to identify all connected components
    void run(const Graph& graph);

    // Accessor for total connected components count
    int getTotalComponents() const;

    // Accessor for sizes of each component
    const std::vector<int>& getComponentSizes() const;

    // Displays formatted BFS component output
    void displayResults() const;
};

#endif // BFS_H
