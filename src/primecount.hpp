#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdint>

#include "phi_block.hpp"

// signum: returns -1, 0, or 1
inline int sgn(const int64_t x)
{
    return (x > 0) - (x < 0);
}


// convenient pi upper bound when x is beyond iacbrtx table
// (Rosser and Schoenfeld 1962)
inline double pi_bound(const uint64_t x)
{
    if (x <= 1) return 1;
    const auto xd = static_cast<double>(x);
    return 1.25506 * xd / log(xd);
}


class Primecount
{
public:
    // tuning parameters
    uint64_t ALPHA;    // tuning parameter (integer here)
    uint64_t C = 8;    // precompute phi_c parameter

    // "global constants"
    uint64_t X;         // Main value to compute pi(X) for
    uint64_t Q;         // phi_c table size, shouldn't be too large
    uint64_t Z;         // X^(2/3) / alpha (approx)
    uint64_t ISQRTX;    // floor(sqrt(X))
    uint64_t IACBRTX;   // floor(alpha cbrt(X))
    uint64_t A;         // pi(alpha cbrt(X))
    uint64_t ASTAR;     // p_a*^2 = alpha cbrt(X), a* = pi(sqrt(alpha cbrt X))
    uint64_t BLOCKMIN;  // minimum block size (in bits)
    uint64_t BLOCKMAX;  // maximum block size (in bits)
    uint64_t K;         // max k for blocks
    uint64_t KL = 64;   // how many k to parallelize at once


    // precomputed tables
    std::vector<int64_t>    MU_PMIN;     // mu(n) pmin(n) for [2,acbrtx]
    std::vector<uint64_t>   PRIMES;      // primes <= acbrtx
    std::vector<uint64_t>   PRIME_COUNT; // pi(x) over [1,acbrtx]
    std::vector<uint64_t>   PHI_C;       // phi(x,c) over [1,Q]
    std::vector<bool>       F_C;         // f(x,c) = [pmin(x) > p_c]
    std::vector<uint64_t>   zks;         // z_k endpoints for interval [1,z]

    Primecount(uint64_t x, uint64_t alpha, uint64_t blockmin, uint64_t blockmax);

    // precompute PRIMES, PRIME_COUNT, MU_PMIN with a standard sieve to acbrtx
    void sieve_mu_prime(size_t SIEVE_SIZE);

    // Precompute phi(x,c), init PHI_C and F_C
    void pre_phi_c();

    // phi(y,c) can be found quickly from the table
    // y can be 10^19 so needs to be unsigned
    uint64_t phi_yc(const uint64_t y) const
    {
        return (y / Q) * PHI_C[Q] + PHI_C[y % Q];
    }

    // contribution of ordinary leaves to phi(x,a)
    uint64_t S0_iter() const;

    // Algorithm 1
    uint64_t S1_iter(uint64_t b, const PhiBlock& phi_block, uint64_t& phi_defer) const;

    // Algorithm 2, reworked from leaves formulas
    uint64_t S2_iter(uint64_t b, const PhiBlock& phi_block, uint64_t& phi_defer) const;

    // Algorithm 3: computation of phi2(x,a)
    uint64_t P2_iter(const PhiBlock& phi_block, uint64_t& v, uint64_t& phi_defer) const;

    // Top-level computation
    uint64_t primecount();
};

