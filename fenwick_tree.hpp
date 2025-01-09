// Credit: cp-algorithms (Jakob Kobler), e-maxx.ru (Maxim Ivanov)

#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>

struct fenwick_tree
{
    size_t len; // 0-based len
    std::vector<int64_t> t; // 1-based tree, indexes [1:len]

    fenwick_tree(const std::vector<bool> &a) :
    len(a.size()) 
    {
        reset(a);
    }

    void reset(const std::vector<bool> &a)
    {
        t.assign(len + 1, 0);
        // linear time construction
        for (size_t i = 1; i <= len; ++i)
        {
            t[i] += a[i-1];
            size_t r = i + (i & -i);
            if (r <= len) t[r] += t[i];
        }
    }

    // 0-based input
    int64_t sum_to(size_t r) const
    {
        int64_t s = 0;
        for (++r; r > 0; r -= r & -r)
            s += t[r];
        return s;
    }

    // 0-based input
    void add_to(size_t i, int64_t delta)
    {
        for (++i; i <= len; i += i & -i)
            t[i] += delta;
    }
};
