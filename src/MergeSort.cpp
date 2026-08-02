#include "MergeSort.h"

using namespace std;

void MergeSort::merge(vector<Gene>& genes, int left, int mid, int right) {
    int n1 = mid - left + 1;
    int n2 = right - mid;

    vector<Gene> L(n1);
    vector<Gene> R(n2);

    for (int i = 0; i < n1; ++i) {
        L[i] = genes[left + i];
    }
    for (int j = 0; j < n2; ++j) {
        R[j] = genes[mid + 1 + j];
    }

    int i = 0;
    int j = 0;
    int k = left;

    while (i < n1 && j < n2) {
        if (L[i].getVariance() >= R[j].getVariance()) {
            genes[k] = L[i];
            i++;
        } else {
            genes[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        genes[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        genes[k] = R[j];
        j++;
        k++;
    }
}

void MergeSort::mergeSortHelper(vector<Gene>& genes, int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;

        mergeSortHelper(genes, left, mid);
        mergeSortHelper(genes, mid + 1, right);

        merge(genes, left, mid, right);
    }
}

void MergeSort::sort(vector<Gene>& genes) {
    if (genes.empty()) return;
    mergeSortHelper(genes, 0, static_cast<int>(genes.size()) - 1);
}
