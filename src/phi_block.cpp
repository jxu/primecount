#include <cassert>
#include "phi_block.hpp"

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// Fits exactly in a power of 2 space!
// Using the MSB for the underlying bool array isn't faster but it's bit fun

FenwickTree::FenwickTree(const std::vector<bool>& ind) : t(ind.size(), 0)
{
    assert(ind.size() < MSB_MASK);

    // fancy linear time construction
    for (uint32_t i = 0; i < t.size(); ++i)
    {
        t[i] += ind[i]; // add input 0/1 to sums
        // should stay with MSB unset
        uint32_t r = i | (i + 1);
        if (r < t.size())
            t[r] += t[i]; // push forward (ignoring MSB)

        t[i] |= ind[i] << 31; // set MSB with input bool
    }
}

// sum values a[0..r] (0-based)
uint32_t FenwickTree::sum_to(uint32_t r) const
{
    assert(r < t.size());

    uint32_t s = 0;
    // r can go "negative" here in the weird 0-indexed tree
    for (; r < MSB_MASK; r = (r & (r + 1)) - 1)
        s += t[r];
    return s & ~MSB_MASK; // return without MSB
}

// will only decrease if underlying ind[i] = 1
// 0-based input
void FenwickTree::try_decrease(uint32_t i)
{
    assert(i < t.size());

    if (t[i] & MSB_MASK) // if set
    {
        t[i] &= ~MSB_MASK; // unset
        for (; i < t.size(); i |= (i + 1))
            --t[i]; // ignore MSB
    }
}


// In the whole computation, Bk processed sequentially for k = 1 to K
// (potentially parallelizable by tracking phi(zk1-1,b) base values used)
// Within each Bk:
// For b from 1 to a:
//   Sieve out pb
//   Update S1b and S2b

// The physical block size is half the logical size by only storing odd values
// For example, [51, 101) would map to ind [0, 25) via y -> (y-zk1)/2

void PhiBlock::sieve_out(uint64_t pb)
{
    assert(pb > 2); // 2 already sieved by default

    // sieve out pb
    if (zk1 <= pb && pb < zk)
        phi_sum.try_decrease(tree_index(pb));

    // now only need to start at pb^2
    // (doesn't really help)
    uint64_t j0 = std::max(pb*pb, pb * ceil_div(zk1, pb));
    if (j0 % 2 == 0)
        j0 += pb; // ensure odd

    for (uint64_t j = j0; j < zk; j += 2*pb)
    {
        phi_sum.try_decrease(tree_index(j));
    }
}