#include "fenwick_tree.hpp"
#include <cassert>

void fenwick_tree::reset()
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

int32_t fenwick_tree::sum_to(size_t r) const
{
    assert(r+1 < t.size());
    // special fenwick tree sum operation
    int32_t s = 0;
    for (++r; r > 0; r -= r & -r)
        s += t[r];
    return s;
}

void fenwick_tree::try_decrease(size_t i)
{   
    assert(i < ind.size());
    if (ind[i])
    {
        ind[i] = 0;
        for (++i; i <= len; i += i & -i)
            --t[i];
    }
}

