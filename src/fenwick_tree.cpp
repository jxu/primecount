#include "fenwick_tree.hpp"

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// Fits exactly in a power of 2 space!
// Using the MSB for the underlying bool array isn't faster but it's bit fun

fenwick_tree::fenwick_tree(const std::vector<bool> &ind) : t(ind.size(), 0)
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
uint32_t fenwick_tree::sum_to(uint32_t r) const
{
    assert(size_t(r) < t.size());

    uint32_t s = 0;
    // r can go "negative" here in the weird 0-indexed tree
    for (; r < MSB_MASK; r = (r & (r + 1)) - 1)
        s += t[r];
    return s & ~MSB_MASK; // return without MSB
}

// will only decrease if underlying ind[i] = 1
// 0-based input
void fenwick_tree::try_decrease(uint32_t i)
{
    assert(i < t.size());

    if (t[i] & MSB_MASK) // if set
    {
        t[i] &= ~MSB_MASK; // unset
        for (; i < t.size(); i |= (i + 1))
            --t[i]; // ignore MSB
    }
}
