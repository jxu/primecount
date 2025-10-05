#pragma once
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdint>

#include "fenwick_tree.hpp"

// ceil(x/y) for positive integers
inline int64_t ceil_div(int64_t x, int64_t y)
{
    return x / y + (x % y > 0);
}

// Represents Bk = [zk1, zk) that partition [1, ceil(z)]
// The interval size should be O(iacbrtx) in theory
//
// In the whole computation, Bk processed sequentially for k = 1 to K
// (potentially parallelizable by tracking phi(zk1-1,b) base values used)
// Within each Bk:
// For b from 1 to a:
//   Sieve out pb
//   Update S1b and S2b

// The physical block size is half the logical size by only storing odd values
// For example, [51, 101) would map to ind [0, 25) via y -> (y-zk1)/2

class PhiBlock
{
public:
    int64_t          zk1;       // z_{k-1}, block lower bound (inclusive)
    int64_t          zk;        // z_k, block upper bound (exclusive)
    size_t           bsize;     // logical block size
    size_t           psize;     // physical block size

    fenwick_tree     phi_sum;   // data structure for efficient partial sums

    // construct new block without k or b explicitly
    PhiBlock(const std::vector<bool>& ind, size_t zk1, size_t zk) :
        zk1(zk1),
        zk(zk),
        bsize(zk - zk1),
        psize(bsize / 2),
        phi_sum(ind)
    {
        assert(zk1 % 2 == 1);
        assert(bsize % 2 == 0);
    }

    // translate into actual index into tree
    inline size_t tree_index(const int64_t y) const
    {
        return (y - zk1)/2;
    }

    // sum contribution of this block to phi(y,b)
    // = phi(y,b) - phi(zk1 - 1, b) base
    int64_t sum_to(int64_t y) const
    {
        assert(y >= zk1);
        assert(y < zk);
        return phi_sum.sum_to(tree_index(y));
    }

    // sieve out multiples of p_b for this block (including p_b)
    void sieve_out(int64_t pb)
    {
        assert(pb > 2); // 2 already sieved by default

        // sieve out pb
        if (zk1 <= pb && pb < zk)
            phi_sum.try_decrease(tree_index(pb));

        // now only need to start at pb^2
        // (doesn't really help)
        int64_t j0 = std::max(pb*pb, pb * ceil_div(zk1, pb));
        if (j0 % 2 == 0)
            j0 += pb; // ensure odd

        for (int64_t j = j0; j < zk; j += 2*pb)
        {
            phi_sum.try_decrease(tree_index(j));
        }
    }
};
