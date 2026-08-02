#include "Queue.h"

using namespace std;

Queue::Queue() : head(0) {}

void Queue::enqueue(int val) {
    data.push_back(val);
}

void Queue::dequeue() {
    if (!isEmpty()) {
        head++;
    }
}

int Queue::front() const {
    if (!isEmpty()) {
        return data[head];
    }
    return -1;
}

bool Queue::isEmpty() const {
    return head >= data.size();
}

size_t Queue::size() const {
    if (isEmpty()) return 0;
    return data.size() - head;
}

void Queue::clear() {
    data.clear();
    head = 0;
}
