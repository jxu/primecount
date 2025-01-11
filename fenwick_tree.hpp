// Credit: cp-algorithms (Jakob Kobler), e-maxx.ru (Maxim Ivanov)

#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>

// customized to save memory 
// only holds unsigned 32-bit values, takes in bools as input
struct fenwick_tree
{
    uint32_t len; // 0-based len
    std::vector<uint32_t> t; // 1-based tree, indexes [1:len]

    // input bool vector
    fenwick_tree(const std::vector<bool> &a) :
    len(a.size()) 
    {
        reset(a);
    }

    void reset(const std::vector<bool> &a)
    {
        t.assign(len + 1, 0);
        // linear time construction
        for (uint32_t i = 1; i <= len; ++i)
        {
            t[i] += a[i-1];
            uint32_t r = i + (i & -i);
            if (r <= len) t[r] += t[i];
        }
    }

    // 0-based input
    uint32_t sum_to(uint32_t r) const
    {
        uint32_t s = 0;
        for (++r; r > 0; r -= r & -r)
            s += t[r];
        return s;
    }

    // 0-based input
    void decrease(uint32_t i)
    {
        for (++i; i <= len; i += i & -i)
            --t[i];
    }
};
