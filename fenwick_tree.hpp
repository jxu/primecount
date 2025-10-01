#pragma once
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdint>

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
class fenwick_tree
{
private:
    size_t                  len; // 0-based len
    std::vector<bool>       ind; // underlying array
    std::vector<uint32_t>   t;   // 1-based tree, indexes [1:len]

public:
    // init with input bool array
    fenwick_tree(const std::vector<bool>& ind) :
        len(ind.size()),
        ind(ind),
        t(len + 1, 0)
    {
        // linear time construction
        for (size_t i = 1; i <= len; ++i)
        {
            t[i] += ind[i-1];
            size_t r = i + (i & -i);
            if (r <= len) 
                t[r] += t[i];
        }
    }

    // sum values a[0..r] (0-based)
    uint32_t sum_to(size_t r) const
    {
        assert(r+1 < t.size());
        uint32_t s = 0;
        for (++r; r > 0; r -= r & -r)
            s += t[r];
        return s;
    }

    // will only decrease if ind[i] = 1
    // 0-based input
    void try_decrease(size_t i)
    {
        assert(i < t.size());
        if (ind[i])
        {
            ind[i] = 0;
            for (++i; i <= len; i += i & -i)
                --t[i];
        }
    }
};

