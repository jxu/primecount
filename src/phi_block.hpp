#pragma once

#include <cassert>
#include "fenwick_tree.hpp"

// ceil(x/y) for positive integers
inline uint64_t ceil_div(uint64_t x, uint64_t y)
{
    return x / y + (x % y > 0);
}


// Represents Bk = [zk1, zk) that partition [1, ceil(z)]
// The interval size should be O(iacbrtx) in theory
class PhiBlock
{
public:
    uint64_t    zk1;       // z_{k-1}, block lower bound (inclusive)
    uint64_t    zk;        // z_k, block upper bound (exclusive)
    size_t      bsize;     // logical block size
    size_t      psize;     // physical block size

    FenwickTree    phi_sum;   // data structure for efficient partial sums

    // construct new block without k or b explicitly
    PhiBlock(const std::vector<bool>& ind, uint64_t zk1, uint64_t zk) :
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
    uint64_t tree_index(const uint64_t y) const
    {
        return (y - zk1)/2;
    }

    // sum contribution of this block to phi(y,b)
    // = phi(y,b) - phi(zk1 - 1, b) base
    uint64_t sum_to(uint64_t y) const
    {
        assert(y >= zk1);
        assert(y < zk);
        return phi_sum.sum_to(tree_index(y));
    }

    // sieve out multiples of p_b for this block (including p_b)
    void sieve_out(uint64_t pb);
};