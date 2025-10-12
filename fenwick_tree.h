#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

const uint32_t MSB_MASK = 1 << 31;

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// Fits exactly in a power of 2 space!
// Using the MSB for the underlying bool array isn't faster but it's bit fun
typedef struct
{
    uint32_t  len;
    uint32_t* t; // 0-indexed tree
} fenwick_tree;

// init with input bool array
fenwick_tree ft_new(bool ind[], uint32_t len)
{
    assert(len < MSB_MASK);
    // init array with zeros
    uint32_t* t = calloc(len, sizeof *t); 

    // fancy linear time construction
    for (uint32_t i = 0; i < len; ++i)
    {
        t[i] += ind[i]; // add input 0/1 to sums
        // should stay with MSB unset
        uint32_t r = i | (i + 1);
        if (r < len)
            t[r] += t[i]; // push forward (ignoring MSB)
        
        t[i] |= ind[i] << 31; // set MSB with input bool
    }

    fenwick_tree ft;
    ft.len = len;
    ft.t = t;

    return ft;
}

// sum values a[0..r] (0-based)
uint32_t ft_sum_to(fenwick_tree ft, int32_t r)
{
    assert(r < ft.len);

    uint32_t s = 0;
    // r can go negative here in the weird 0-indexed tree
    for (; r >= 0; r = (r & (r + 1)) - 1)
        s += ft.t[r];
    return s & ~MSB_MASK; // return without MSB
}

// will only decrease if underlying ind[i] = 1
// 0-based input
// ft is passed by value, but underlying t is the same
void ft_try_decrease(fenwick_tree ft, uint32_t i)
{
    assert (i < ft.len);

    if (ft.t[i] & MSB_MASK) // if set
    {
        ft.t[i] &= ~MSB_MASK; // unset
        for (; i < ft.len; i |= (i + 1))
            --ft.t[i]; // ignore MSB
    }
}

