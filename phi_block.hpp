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
        assert(ind.size() == psize);
        assert(bsize % 2 == 0);
    }

    // translate into actual index into tree
    size_t tree_index(const int64_t y) const
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
    void sieve_out(size_t pb)
    {
        assert(pb > 2); // 2 already sieved by default

        int64_t jstart = pb * ceil_div(zk1, pb);
        if (jstart % 2 == 0)
            jstart += pb; // ensure odd

        for (int64_t j = jstart; j < zk; j += 2*pb)
        {
            phi_sum.try_decrease(tree_index(j));
        }
    }
};
