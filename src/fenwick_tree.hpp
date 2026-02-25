#pragma once
#include <vector>
#include <cassert>
#include <cstdint>


class FenwickTree
{
private:
    const uint32_t MSB_MASK = 1UL << 31;
    std::vector<uint32_t>   t;   // 0-indexed tree

public:
    // init with input bool array
    FenwickTree(const std::vector<bool>& ind);

    // sum values a[0..r] (0-based)
    uint32_t sum_to(uint32_t r) const;

    // will only decrease if underlying ind[i] = 1
    // 0-based input
    void try_decrease(uint32_t i);
};
