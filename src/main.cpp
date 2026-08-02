#include "../include/ProjectManager.h"
#include <cstdlib>
#include <ctime>

using namespace std;

// Program entry point
int main() {
    // Seed RNG once at program startup so every run produces different results
    srand(static_cast<unsigned int>(time(NULL)));

    // Instantiate project manager controller
    ProjectManager manager;
    
    // Start main application loop
    manager.start();
    
    return 0;
}
