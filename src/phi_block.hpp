#pragma once

#include <vector>
#include "fenwick_tree.hpp"

using namespace std; // bad in header

// ceil(x/y) for positive integers
inline long ceil_div(long x, long y)
{
    return x / y + (x % y > 0);
}


// Represents Bk = [zk1, zk) and functions to compute phi(y,b)
// phi(y,b) = count (pos) integers <= y not divisible by the first b primes
// i.e. not p_b-smooth
// The physical block size is half the logical size by only storing odd values
// For example, [51, 101) would map to ind [0, 25) via y -> (y-zk1)/2

struct PhiBlock
{
    size_t           bsize;     // logical block size
    size_t           psize;     // physical block size
    size_t           zk1;       // z_{k-1}, block lower bound (inclusive)
    size_t           zk;        // z_k, block upper bound (exclusive)
    fenwick_tree     phi_sum;   // data structure for efficient partial sums
    
    // phi_save(k,b) = phi(zk1-1,b) from prev block, for implicit b
    // c <= b < a-1
    vector<long>     phi_save;  

    // init block at k=1
    PhiBlock(size_t a, size_t bsize) :
        bsize(bsize),
        psize(bsize / 2),  // actually only store odd values
        phi_sum(psize),
        phi_save(a + 1, 0)
    {
        assert(bsize % 2 == 0); // only tested even sizes
    }

    void new_block(size_t k)
    {
        // TODO: more flexible zk1 and zk?
        zk1 = bsize * (k-1) + 1;
        zk = bsize * k + 1;
        phi_sum.reset();
        // does not reset phi_save!
        
        //cout << "block [" << zk1 << "," << zk << ")\n";
    }

    // translate into actual index into tree
    size_t tree_index(const size_t y) const
    {
        return (y - zk1)/2;
    }

    // phi(y,b) compute 
    long sum_to(size_t y, size_t b) const
    {
        assert(y >= zk1);
        assert(y < zk);
        return phi_save[b] + phi_sum.sum_to(tree_index(y));
    }

    // sieve out p_b for this block
    void sieve_out(size_t pb)
    {
        assert(pb > 2);

        size_t jstart = pb * ceil_div(zk1, pb);
        if (jstart % 2 == 0) 
            jstart += pb; // ensure odd
        
        for (size_t j = jstart; j < zk; j += 2*pb)
        {
            phi_sum.try_decrease(tree_index(j)); 
        }
    }

    void update_save(size_t b)
    {
        phi_save[b] += phi_sum.sum_to(psize-1);
    }
};


