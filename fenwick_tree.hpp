#pragma once
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdint>
#include <stdexcept>

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// fits exactly in a power of 2 space!
class fenwick_tree
{
private:
    std::vector<bool>       ind; // underlying array
    std::vector<uint32_t>   t;   // 0-indexed tree

public:
    // init with input bool array
    fenwick_tree(const std::vector<bool>& ind) :
        ind(ind),
        t(ind.size(), 0)
    {
        // linear time construction
        for (size_t i = 0; i < t.size(); ++i)
        {
            t[i] += ind[i];
            size_t r = i | (i + 1);
            if (r < t.size()) 
                t[r] += t[i];
        }
    }

    // sum values a[0..r] (0-based)
    uint32_t sum_to(int32_t r) const
    {
        if ((size_t)r >= t.size())
            throw std::out_of_range("sum_to out of range");

        uint32_t s = 0;
        for (; r >= 0; r = (r & (r + 1)) - 1)
            s += t[r];
        return s;
    }

    // will only decrease if ind[i] = 1
    // 0-based input
    void try_decrease(size_t i)
    {
        if (i >= t.size()) 
            throw std::out_of_range("try_decrease out of range");

        if (ind[i])
        {
            ind[i] = 0;
            for (; i < t.size(); i |= (i + 1))
                --t[i];
        }
    }
};

