#ifndef QUEUE_H
#define QUEUE_H

#include <vector>
#include <cstddef>

// Custom FIFO Queue implementation using std::vector container
class Queue {
private:
    std::vector<int> data;
    size_t head;

public:
    Queue();

    // Adds element to the back of the queue
    void enqueue(int val);

    // Removes element from the front of the queue
    void dequeue();

    // Returns front element without removing it (-1 if empty)
    int front() const;

    // Checks if queue is empty
    bool isEmpty() const;

    // Returns number of elements in queue
    size_t size() const;

    // Resets queue
    void clear();
};

#endif // QUEUE_H
