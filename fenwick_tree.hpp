#pragma once

#include <vector>
#include <cstddef>


struct fenwick_tree {
    size_t len; // 0-based len
    std::vector<long long> t; // 1-based tree, indexes [1:len]

    fenwick_tree(std::vector<long long> const &a) 
    {
        len = a.size();
        t.assign(len + 1, 0);
        for (size_t i = 0; i < len; ++i)
            add_to(i, a[i]);
    }

    long long sum_to(size_t r) 
    {
        long long s = 0;
        for (++r; r > 0; r -= r & -r)
            s += t[r];
        return s;
    }

    void add_to(size_t i, long long delta) 
    {
        for (++i; i <= len; i += i & -i)
            t[i] += delta;   
    }
};
