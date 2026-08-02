#ifndef STACK_H
#define STACK_H

#include <vector>
#include <cstddef>

// Custom LIFO Stack implementation using std::vector container
class Stack {
private:
    std::vector<int> data;

public:
    Stack();

    // Pushes element to the top of the stack
    void push(int val);

    // Removes top element from the stack
    void pop();

    // Returns top element without removing it (-1 if empty)
    int top() const;

    // Checks if stack is empty
    bool isEmpty() const;

    // Returns number of elements in stack
    size_t size() const;

    // Resets stack
    void clear();
};

#endif // STACK_H
