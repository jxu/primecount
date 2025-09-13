#pragma once

#include <cstdint>
#include <vector>
#include <cassert>

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
    void reset()
    {
        ind.assign(len + 1, 1);
        t.assign(len + 1, 0);
        // fancy linear time construction
        for (size_t i = 1; i <= len; ++i)
        {
            t[i] += ind[i-1];
            size_t r = i + (i & -i);
            if (r <= len) t[r] += t[i];
        }
    }
    
    // sum values a[0..r] (0-based)
    int32_t sum_to(size_t r) const
    {
        assert(r+1 < t.size());
        // special fenwick tree sum operation
        int32_t s = 0;
        for (++r; r > 0; r -= r & -r)
            s += t[r];
        return s;
    }

    // will only decrease if ind[i] is not already marked
    // 0-based input
    void try_decrease(size_t i)
    {   
        assert(i < ind.size());
        if (ind[i])
        {
            ind[i] = 0;
            for (++i; i <= len; i += i & -i)
                --t[i];
        }
    }
};

