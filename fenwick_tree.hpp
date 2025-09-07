#pragma once

#include <cstdint>
#include <vector>

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// only holds 32-bit values
class fenwick_tree
{
private:
    size_t                  len; // 0-based len
    std::vector<int32_t>    t;   // 1-based tree, indexes [1:len]
    std::vector<bool>       ind; // 0/1 to track [pmin(y) > pb]

public:
    // init array of 1s of length psize
    fenwick_tree(size_t psize) : len(psize)
    {
        reset();
    }

    // reset ind to all 1s and reconstruct t
    void reset();

    // sum values a[0..r] (0-based)
    int32_t sum_to(size_t r) const;

    // will only decrease if ind[i] is not already marked
    // 0-based input
    void try_decrease(size_t i);
};


