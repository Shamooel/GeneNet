#ifndef MERGESORT_H
#define MERGESORT_H

#include <vector>
#include "Gene.h"

// Class implementing manual Merge Sort algorithm for sorting Genes by variance
class MergeSort {
private:
    // Helper function to merge two sorted subarrays in descending order of variance
    static void merge(std::vector<Gene>& genes, int left, int mid, int right);

    // Recursive helper function for divide-and-conquer merge sort
    static void mergeSortHelper(std::vector<Gene>& genes, int left, int right);

public:
    // Sorts vector of Genes by variance in descending order (highest variance first)
    static void sort(std::vector<Gene>& genes);
};

#endif // MERGESORT_H
