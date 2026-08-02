#include "Menu.h"
#include <iostream>

using namespace std;

void Menu::displayMainMenu() const {
    cout << "========================================\n";
    cout << "GeneNet Project\n";
    cout << "========================================\n";
    cout << "1. Run Complete Pipeline\n";
    cout << "2. Run Individual Module\n";
    cout << "3. Delete Results\n";
    cout << "4. Exit\n";
    cout << "Enter choice: ";
}

void Menu::displayIndividualModuleMenu() const {
    cout << "========================================\n";
    cout << "Individual Modules\n";
    cout << "========================================\n";
    cout << "1. CSV Reader\n";
    cout << "2. Variance\n";
    cout << "3. Top Gene Selection\n";
    cout << "4. Hash Table\n";
    cout << "5. Pearson Correlation\n";
    cout << "6. Graph\n";
    cout << "7. BFS\n";
    cout << "8. DFS\n";
    cout << "9. Decision Tree\n";
    cout << "10. Bagging\n";
    cout << "11. Evaluation\n";
    cout << "0. Back\n";
    cout << "Enter choice: ";
}

void Menu::displayDeleteMenu() const {
    cout << "========================================\n";
    cout << "Delete Results\n";
    cout << "========================================\n";
    cout << "1. Delete Complete Pipeline Results\n";
    cout << "2. Delete Individual Module Results\n";
    cout << "3. Delete Everything\n";
    cout << "0. Back\n";
    cout << "Enter choice: ";
}

int Menu::getChoice() const {
    int choice = 0;
    if (!(cin >> choice)) {
        cin.clear();
        cin.ignore(10000, '\n');
        return -1;
    }
    return choice;
}

bool Menu::askViewResults() const {
    cout << "\nDo you want to view the results?\n";
    cout << "1. Yes\n";
    cout << "2. No\n";
    cout << "Enter choice: ";
    int c = getChoice();
    return (c == 1);
}

void Menu::displayMessage(const string& message) const {
    cout << message << "\n";
}
