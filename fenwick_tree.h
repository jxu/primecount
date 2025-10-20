#pragma once
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

typedef uint32_t* fenwick_tree; // light abstraction

const uint32_t MSB_MASK = 1 << 31;

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// Using the MSB for the underlying bool array isn't faster but it's bit fun

// init with input bool array
// result is 1-indexed int array of size len+1, 0-index stores size
fenwick_tree ft_new(bool* ind, uint32_t len)
{
    assert(len + 1 < MSB_MASK);
    // init array with zeros
    fenwick_tree ft = calloc(len + 1, sizeof *ft);

    ft[0] = len;

    // fancy linear time construction
    for (uint32_t i = 1; i <= len; ++i)
    {
        ft[i] += ind[i-1]; // add input 0/1 to sums
        // should stay with MSB unset
        uint32_t r = i + (i & -i);
        if (r <= len)
            ft[r] += ft[i]; // push forward (ignoring MSB)
        
        ft[i] |= ind[i-1] << 31; // set MSB with input bool
    }

    return ft;
}

// sum values a[0..r] (0-based)
uint32_t ft_sum_to(fenwick_tree ft, uint32_t r)
{
    assert(r < ft[0]); // TODO: make into function

    uint32_t s = 0;
    // r can go negative here in the weird 0-indexed tree
    for (++r; r > 0; r -= r & -r)
        s += ft[r];
    return s & ~MSB_MASK; // return without MSB
}

// will only decrease if underlying ind[i] = 1
// 0-based input
// ft is passed by value, but underlying t is the same
void ft_try_decrease(fenwick_tree ft, uint32_t i)
{
    assert (i < ft[0]);
    ++i;

    if (ft[i] & MSB_MASK) // if set
    {
        ft[i] &= ~MSB_MASK; // unset
        for (; i <= ft[0]; i += i & -i)
            --ft[i]; // ignore MSB
    }
}

void ft_delete(fenwick_tree ft)
{
    free(ft);
}
