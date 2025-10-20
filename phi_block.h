#pragma once
#include <assert.h>
#include <stdint.h>

#include "fenwick_tree.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))


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
//   For b from 1 to a:
//     Sieve out pb
//     Update S1b and S2b

// The physical block size is half the logical size by only storing odd values
// For example, [51, 101) would map to ind [0, 25) via y -> (y-zk1)/2

typedef struct
{
    int64_t          zk1;       // z_{k-1}, block lower bound (inclusive)
    int64_t          zk;        // z_k, block upper bound (exclusive)
    size_t           bsize;     // logical block size
    size_t           psize;     // physical block size

    fenwick_tree     phi_sum;   // data structure for efficient partial sums
} phi_block;


// construct new block without k or b explicitly
// ind should be length zk - zk1
phi_block* phi_block_new(bool ind[], size_t zk1, size_t zk)
{
    phi_block* block = malloc(sizeof(phi_block));

    block->zk1 = zk1;
    block->zk = zk;
    block->bsize = zk - zk1;
    block->psize = block->bsize / 2;

    assert(zk1 % 2 == 1);
    assert(block->bsize % 2 == 0);

    block->phi_sum = ft_new(ind, block->bsize);

    return block;
}

// translate into actual index into tree
inline size_t tree_index(phi_block* block, const int64_t y)
{
    return (y - block->zk1)/2;
}

// sum contribution of this block to phi(y,b)
// = phi(y,b) - phi(zk1 - 1, b) base
int64_t sum_to(phi_block* block, int64_t y)
{
    assert(y >= block->zk1);
    assert(y < block->zk);
    return ft_sum_to(block->phi_sum, tree_index(block, y));
}

// sieve out multiples of p_b for this block (including p_b)
void sieve_out(phi_block* block, int64_t pb)
{
    assert(pb > 2); // 2 already sieved by default

    // sieve out pb
    if (block->zk1 <= pb && pb < block->zk)
        ft_try_decrease(block->phi_sum, tree_index(block, pb));

    // now only need to start at pb^2
    // (doesn't really help)
    int64_t j0 = MAX(pb*pb, pb * ceil_div(block->zk1, pb));
    if (j0 % 2 == 0)
        j0 += pb; // ensure odd

    for (int64_t j = j0; j < block->zk; j += 2*pb)
    {
        ft_try_decrease(block->phi_sum, tree_index(block, j));
    }
}

void phi_block_delete(phi_block* block)
{
    ft_delete(block->phi_sum);
    free(block);
}
