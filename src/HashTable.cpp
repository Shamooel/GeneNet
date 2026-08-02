#include "HashTable.h"

using namespace std;

HashTable::HashTable(size_t cap) : capacity(cap), size(0) {
    table.resize(capacity);
}

size_t HashTable::hashFunction(const string& key) const {
    unsigned long hash = 5381;
    for (size_t i = 0; i < key.length(); ++i) {
        hash = ((hash << 5) + hash) + static_cast<unsigned char>(key[i]);
    }
    return hash % capacity;
}

void HashTable::insert(const string& key, int value) {
    size_t index = hashFunction(key);
    for (size_t i = 0; i < table[index].size(); ++i) {
        if (table[index][i].key == key) {
            table[index][i].value = value;
            return;
        }
    }
    table[index].push_back(HashNode(key, value));
    size++;
}

int HashTable::get(const string& key) const {
    size_t index = hashFunction(key);
    for (size_t i = 0; i < table[index].size(); ++i) {
        if (table[index][i].key == key) {
            return table[index][i].value;
        }
    }
    return -1;
}

bool HashTable::contains(const string& key) const {
    return get(key) != -1;
}

void HashTable::clear() {
    for (size_t i = 0; i < table.size(); ++i) {
        table[i].clear();
    }
    size = 0;
}

size_t HashTable::getSize() const {
    return size;
}

size_t HashTable::getCapacity() const {
    return capacity;
}
