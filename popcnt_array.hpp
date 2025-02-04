#pragma once

#include <vector>
#include <cstdint>

#include <cassert>
#include <bit>


struct popcnt_array
{
    uint32_t len; // length of array in 64-bit words
    std::vector<uint64_t> ind;  //  0/1 to track [pmin(y) > pb]
    std::vector<uint32_t> t; // t[i] = sums array[0:i/64)

    // for now, l is limited to 32-bit size so the total fits in 32-bit
    popcnt_array(const uint32_t l) : 
        len(l / 64 + 1), // ith index belongs to word w_(i/64), round up
        ind(len, ~0), // all 1s
        t(len)
    {
        reset();       
    }

    void reset()
    {
        ind.assign(len, ~0);
        
        refresh();
    }
    
    // sum array[0:y) exclusive!
    uint64_t sum_to(const uint32_t y) const 
    {
        uint64_t mask = (1ULL << (y % 64)) - 1;
        return t[y/64] + std::popcount(ind[y/64] & mask);
    }

    // unset index i 
    void unset(const uint32_t y)
    {
        uint64_t mask = (1ULL << (y % 64)); // (y % 64)th bit set
        ind[y / 64] &= ~mask;
        // does not update t 
    }

    // sync t with ind
    void refresh(void)
    {
        t[0] = 0;
        for (uint32_t i = 1; i < len; ++i)
        {
            t[i] = std::popcount(ind[i-1]) + t[i-1];
        }
    }

};

