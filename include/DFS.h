#ifndef DFS_H
#define DFS_H

#include <vector>
#include "Graph.h"
#include "Stack.h"

// Class for Iterative Depth-First Search connected component traversal using custom Stack
class DFS {
private:
    std::vector<int> componentSizes;

public:
    DFS();

    // Runs iterative DFS traversal on graph to identify all connected components
    void run(const Graph& graph);

    // Accessor for total connected components count
    int getTotalComponents() const;

    // Accessor for sizes of each component
    const std::vector<int>& getComponentSizes() const;

    // Displays formatted DFS component output
    void displayResults() const;
};

#endif // DFS_H
