#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <string>
#include <vector>

// Struct representing a Key-Value pair node in the Hash Table
struct HashNode {
    std::string key;
    int value;

    HashNode() : key(""), value(-1) {}
    HashNode(const std::string& k, int v) : key(k), value(v) {}
};

// Custom Hash Table class using separate chaining with vector buckets
class HashTable {
private:
    std::vector<std::vector<HashNode> > table;
    size_t capacity;
    size_t size;

    // DJB2 polynomial string hash function
    size_t hashFunction(const std::string& key) const;

public:
    HashTable(size_t cap = 1009);

    // Inserts key-value pair into hash table
    void insert(const std::string& key, int value);

    // Retrieves value for key (-1 if not found)
    int get(const std::string& key) const;

    // Checks if key exists in hash table
    bool contains(const std::string& key) const;

    // Resets/clears the hash table
    void clear();

    // Returns total stored entries count
    size_t getSize() const;

    // Returns table capacity
    size_t getCapacity() const;
};

#endif // HASHTABLE_H
