#pragma once
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdint>
#include <stdexcept>

const uint32_t MSB_MASK = 1 << 31;

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// Fits exactly in a power of 2 space!
// Using the MSB for the underlying bool array isn't faster but it's bit fun
class fenwick_tree
{
private:
    std::vector<uint32_t>   t;   // 0-indexed tree

public:
    // init with input bool array
    fenwick_tree(const std::vector<bool>& ind) :
        t(ind.size(), 0)
    {
        if (ind.size() > MSB_MASK)
            throw std::out_of_range("input vec is limited to 31 bits");

        // fancy linear time construction
        for (uint32_t i = 0; i < t.size(); ++i)
        {
            t[i] += ind[i]; // add input 0/1 to sums
            // should stay with MSB unset
            uint32_t r = i | (i + 1);
            if (r < t.size())
                t[r] += t[i]; // push forward (ignoring MSB)
            
            t[i] |= ind[i] << 31; // set MSB with input bool
        }
    }

    // sum values a[0..r] (0-based)
    uint32_t sum_to(int32_t r) const
    {
        if (size_t(r) >= t.size())
            throw std::out_of_range("sum_to out of range");

        uint32_t s = 0;
        // r can go negative here in the weird 0-indexed tree
        for (; r >= 0; r = (r & (r + 1)) - 1)
            s += t[r];
        return s & ~MSB_MASK; // return without MSB
    }

    // will only decrease if underlying ind[i] = 1
    // 0-based input
    void try_decrease(uint32_t i)
    {
        if (i >= t.size())
            throw std::out_of_range("try_decrease out of range");

        if (t[i] & MSB_MASK) // if set
        {
            t[i] &= ~MSB_MASK; // unset
            for (; i < t.size(); i |= (i + 1))
                --t[i]; // ignore MSB
        }
    }
};

