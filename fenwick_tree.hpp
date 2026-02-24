#pragma once
#include <vector>
#include <cassert>
#include <cstdint>

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// Fits exactly in a power of 2 space!
// Using the MSB for the underlying bool array isn't faster but it's bit fun
class fenwick_tree
{
private:
    const uint32_t MSB_MASK = 1UL << 31;
    std::vector<uint32_t>   t;   // 0-indexed tree

public:
    // init with input bool array
    fenwick_tree(const std::vector<bool>& ind);

    // sum values a[0..r] (0-based)
    uint32_t sum_to(uint32_t r) const;

    // will only decrease if underlying ind[i] = 1
    // 0-based input
    void try_decrease(uint32_t i);
};
