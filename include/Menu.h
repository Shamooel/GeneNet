#ifndef MENU_H
#define MENU_H

#include <string>

// Class responsible for rendering CLI menus and handling user input prompts
class Menu {
public:
    // Displays main menu options
    void displayMainMenu() const;

    // Displays individual modules submenu options
    void displayIndividualModuleMenu() const;

    // Displays delete results submenu options
    void displayDeleteMenu() const;

    // Reads user integer choice safely
    int getChoice() const;

    // Asks user if they want to view results (1. Yes, 2. No)
    bool askViewResults() const;

    // Utility to display messages
    void displayMessage(const std::string& message) const;
};

#endif // MENU_H
