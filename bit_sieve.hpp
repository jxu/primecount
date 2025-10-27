#pragma once
#include <vector>
#include <cassert>
#include <cstdint>
#include <stdexcept>

// ceil(x/y) for positive integers
inline uint64_t ceil_div(uint64_t x, uint64_t y)
{
    return x / y + (x % y > 0);
}

// Special bit sieve, optimized for range sums with popcount on 64-bit ints
// Idea stolen from Kim Walisch's primecount
// slower asymptotically because I didn't implement better counters
// This class doesn't handle any of the fancy indexing of phi
// fun bit tricks!
class bit_sieve
{
private:
    size_t len;
    std::vector<uint64_t> t;

public:
    int64_t sum; // track sum of 1s

    // init with input bool array
    bit_sieve(const std::vector<bool>& ind) :
        len(ind.size()),
        t((len - 1) / 64 + 1, 0),
        sum(0)
    {
        for (size_t i = 0; i < ind.size(); ++i)
        {
            if (ind[i])
            {
                ++sum;
                t[i / 64] |= 1ULL << (i % 64);
            }
        }
    }

    // sums range [l,r]
    int64_t sum_range(size_t l, size_t r) const
    {
        assert(l <= r);
        assert(r < len);
        uint64_t lw = l / 64, li = l % 64;
        uint64_t rw = r / 64, ri = r % 64;

        uint64_t lmask = ~0ULL << li; // all 1s except lowest l bits
        uint64_t rmask = ~0ULL >> (63 - ri); // all 0s except lowest r+1 bits
        // WRONG: (1ULL << (ri+1)) - 1 as shifting by 64 is UB!
 
        // special case: within one word
        if (lw == rw)
        {
            return std::popcount(t[lw] & lmask & rmask);
        }

        // spanning multiple words
        int64_t s = std::popcount(t[lw] & lmask);

        for (size_t i = lw + 1; i < rw; ++i)
            s += std::popcount(t[i]);

        s += std::popcount(t[rw] & rmask);
        
        return s;

    }

    // 0-based input
    void try_decrease(uint32_t i)
    {
        if (i >= len)
            throw std::out_of_range("try_decrease out of range");

        // decrease sum if bit is set
        if (t[i / 64] & (1ULL << (i % 64)))
            --sum;

        // clear bit
        t[i / 64] &= ~(1ULL << (i % 64));

    }
};

