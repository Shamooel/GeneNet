#include "Stack.h"

using namespace std;

Stack::Stack() {}

void Stack::push(int val) {
    data.push_back(val);
}

void Stack::pop() {
    if (!isEmpty()) {
        data.pop_back();
    }
}

int Stack::top() const {
    if (!isEmpty()) {
        return data.back();
    }
    return -1;
}

bool Stack::isEmpty() const {
    return data.empty();
}

size_t Stack::size() const {
    return data.size();
}

void Stack::clear() {
    data.clear();
}
