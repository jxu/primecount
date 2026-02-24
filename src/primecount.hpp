#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdint>

#include "fenwick_tree.hpp"



// ceil(x/y) for positive integers
inline uint64_t ceil_div(uint64_t x, uint64_t y)
{
    return x / y + (x % y > 0);
}

// Represents Bk = [zk1, zk) that partition [1, ceil(z)]
// The interval size should be O(iacbrtx) in theory
//
// In the whole computation, Bk processed sequentially for k = 1 to K
// (potentially parallelizable by tracking phi(zk1-1,b) base values used)
// Within each Bk:
// For b from 1 to a:
//   Sieve out pb
//   Update S1b and S2b

// The physical block size is half the logical size by only storing odd values
// For example, [51, 101) would map to ind [0, 25) via y -> (y-zk1)/2

class PhiBlock
{
public:
    uint64_t     zk1;       // z_{k-1}, block lower bound (inclusive)
    uint64_t     zk;        // z_k, block upper bound (exclusive)
    size_t  bsize;     // logical block size
    size_t  psize;     // physical block size

    fenwick_tree    phi_sum;   // data structure for efficient partial sums

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
    void sieve_out(uint64_t pb)
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
};

// signum: returns -1, 0, or 1
inline int sgn(int64_t x)
{
    return (x > 0) - (x < 0);
}


// convenient pi upper bound when x is beyond iacbrtx table
// (Rosser and Schoenfeld 1962)
inline double pi_bound(uint64_t x)
{
    if (x <= 1) return 1;
    return 1.25506 * x / log(x);
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
    std::vector<int64_t>        MU_PMIN;     // mu(n) pmin(n) for [2,acbrtx]
    std::vector<uint64_t>        PRIMES;      // primes <= acbrtx
    std::vector<uint64_t>        PRIME_COUNT; // pi(x) over [1,acbrtx]
    std::vector<uint64_t>        PHI_C;       // phi(x,c) over [1,Q]
    std::vector<bool>       F_C;         // f(x,c) = [pmin(x) > p_c]
    std::vector<uint64_t>        zks;         // z_k endpoints for interval [1,z]

    Primecount(uint64_t x, uint64_t alpha, uint64_t blockmin, uint64_t blockmax) :
        ALPHA(alpha),
        X(x),
        // Q after sieve
        // hope floating point truncated values are exact floors
        Z(cbrt(X) * cbrt(X) / ALPHA), // approx
        ISQRTX(sqrt(X)),
        IACBRTX(ALPHA * cbrt(X)),
        BLOCKMIN(blockmin),
        BLOCKMAX(blockmax)
    {
        if (X < 10)
            throw std::invalid_argument("Program not designed for tiny inputs");

        // check alpha isn't set too large
        if (ALPHA > pow(X, 1/6.))
            throw std::invalid_argument("Alpha set too large");

        std::cout << "Z = " << Z << std::endl;
        std::cout << "IACBRTX = " << IACBRTX << std::endl;
        std::cout << "ISQRTX = " << ISQRTX << std::endl;

        // precompute PRIMES, PRIME_COUNT, MU_PMIN
        // Since p_{a+1} may be needed in S2, leave margin
        size_t SIEVE_SIZE = IACBRTX + 200;
        sieve_mu_prime(SIEVE_SIZE);

        A = PRIME_COUNT[IACBRTX];
        std::cout << "a = " << A << std::endl;

        assert(PRIMES.size() > (size_t)A + 1); // need p_{a+1}

        ASTAR = PRIME_COUNT[uint64_t(sqrt(ALPHA) * pow(X, 1/6.))];
        std::cout << "a* = " << ASTAR << std::endl;

        C = std::min(ASTAR, C);
        std::cout << "C = " << C << std::endl;

        // precompute PHI_C tables
        pre_phi_c();

        // create z_k endpoints
        const size_t BSIZE = 1 << BLOCKMAX;
        zks = {1};

        for (uint64_t i = BLOCKMIN; i < BLOCKMAX; ++i)
            zks.push_back((1 << i) + 1);

        for (uint64_t i = 1 + BSIZE; i <= Z + BSIZE; i += BSIZE)
            zks.push_back(i);

        K = zks.size() - 1;
        std::cout << "K = " << K << std::endl;
    }

    // precompute PRIMES, PRIME_COUNT, MU_PMIN with a standard sieve to acbrtx
    void sieve_mu_prime(size_t SIEVE_SIZE)
    {
        MU_PMIN.assign(SIEVE_SIZE+1, 1);    // init to 1s
        MU_PMIN[1] = 1000;                  // define pmin(1) = +inf
        PRIMES.push_back(1);                // p0 = 1 by convention
        PRIME_COUNT.resize(SIEVE_SIZE+1);   // init values don't matter here

        uint64_t pc = 0; // prime counter

        // sieve of Eratosthenes, modification to keep track of mu sign and pmin
        for (size_t j = 2; j <= SIEVE_SIZE; ++j)
        {
            if (MU_PMIN[j] == 1)   // unmarked, so it is prime
            {
                for (size_t i = j; i <= SIEVE_SIZE; i += j)
                {
                    MU_PMIN[i] = (MU_PMIN[i] == 1) ? -j : -MU_PMIN[i];
                }
            }
        }

        // complete MU_PMIN, compute PRIMES and PRIME_COUNT
        for (int64_t j = 2; j <= int64_t(SIEVE_SIZE); ++j)
        {
            if (MU_PMIN[j] == -j)   // prime
            {
                PRIMES.push_back(uint64_t(j));
                ++pc;

                // mark multiples of p^2 as 0 for mu
                for (size_t i = j*j; i <= SIEVE_SIZE; i += j*j)
                    MU_PMIN[i] = 0;
            }

            PRIME_COUNT[j] = pc;
        }
    }

    // Precompute phi(x,c)
    void pre_phi_c()
    {
        // compute Q as product of first C primes
        Q = 1;
        for (size_t i = 1; i <= C; ++i)
            Q *= PRIMES[i];

        PHI_C.resize(Q+1, 0);
        F_C.resize(Q+1, 1); // index up to Q, inclusive
        F_C[0] = 0;

        for (size_t i = 1; i <= C; ++i)
        {
            uint64_t p = PRIMES[i]; // ith prime, mark multiples as 0
            for (uint64_t j = p; j <= Q; j += p)
                F_C[j] = 0;
        }

        // accumulate
        for (size_t i = 1; i <= Q; ++i)
            PHI_C[i] = F_C[i] + PHI_C[i-1];
    }

    // phi(y,c) can be found quickly from the table
    // y can be 10^19 so needs to be unsigned
    uint64_t phi_yc(uint64_t y)
    {
        return (y / Q) * PHI_C[Q] + PHI_C[y % Q];
    }

    // contribution of ordinary leaves to phi(x,a)
    uint64_t S0_iter(void)
    {
        uint64_t S0 = 0;
        for (size_t n = 1; n <= size_t(IACBRTX); ++n)
        {
            if (uint64_t(std::abs(MU_PMIN[n])) > PRIMES[C])
            {
                S0 += sgn(MU_PMIN[n]) * phi_yc(X / n);
            }
        }

        return S0;
    }

    // Algorithm 1
    uint64_t S1_iter(const uint64_t b, const PhiBlock& phi_block, uint64_t& phi_defer)
    {
        uint64_t S1 = 0;
        uint64_t pb1 = PRIMES[b+1];
        uint64_t defer = 0;
        // m decreasing, modified with fixed upper bound
        uint64_t mb = std::min(IACBRTX, X / (phi_block.zk1 * pb1));

        for (; mb * pb1 > IACBRTX; --mb)
        {
            uint64_t y = X / (mb * pb1);

            assert(y >= phi_block.zk1);
            if (y >= phi_block.zk) break;

            if (uint64_t(std::abs(MU_PMIN[mb])) > pb1)
            {
                S1 -= sgn(MU_PMIN[mb]) * phi_block.sum_to(y);
                defer -= sgn(MU_PMIN[mb]);
            }
        }
        phi_defer = defer;
        return S1;
    }

    // Algorithm 2, reworked from leaves formulas
    uint64_t S2_iter(const uint64_t b, const PhiBlock& phi_block, uint64_t& phi_defer)
    {
        uint64_t S2b = 0;
        uint64_t pb1 = PRIMES[b+1];
        uint64_t zk1 = phi_block.zk1;
        uint64_t zk = phi_block.zk;
        uint64_t defer = 0;

        uint64_t d = A; // fixed starting point

        // attempt to optimize starting d bound
        // pd <= x / zk1
        d = std::min(d, uint64_t(pi_bound(X / (pb1 * zk1))));
        // non-trivial leaves should satisfy pd <= max(x/pb1^2, pb1)
        d = std::min(d, uint64_t(pi_bound(X / (pb1 * pb1))));

        for (; d > b + 1; --d)
        {
            uint64_t pd = PRIMES[d];
            uint64_t y = X / (pb1 * pd);

            // it is possible to avoid integer division until hard leaves,
            // but the comparisons need __int128 for overflow on X=1e19
            if (y < zk1)
                continue;
            if (y >= zk)
                break;

            // trivial leaves, should be skipped since already counted
            if (X / (pb1 * pb1) < pd) // X / pb1^2 < pd
                continue;

            // hard leaves
            if (y >= IACBRTX)
            {
                S2b += phi_block.sum_to(y);
                defer++;
            }

            // easy leaves
            // can't use clustered leaves because of setting d = d' jumps
            else
            {
                // sparse easy leaves
                S2b += PRIME_COUNT[y] - b + 1;
            }
        }

        // minimize updating the references
        phi_defer = defer;
        return S2b;
    }

    // Algorithm 3: computation of phi2(x,a)
    // Modification to original algorithm: use aux sieve [u,w] limited to
    // y's range in [zk1,zk) rather than acbrtx size
    uint64_t P2_iter(const PhiBlock& phi_block, uint64_t& v, uint64_t& phi_defer)
    {
        uint64_t P2 = 0;
        uint64_t defer = 0;
        uint64_t v_defer = 0;

        // maintain aux sieve, u = tracking pb starting at max, y = x / u
        // x/zk < u <= x/zk1
        // iacbrtx < u <= isqrtx
        uint64_t u = std::min(ISQRTX, X / phi_block.zk1);
        uint64_t w = std::max(IACBRTX, X / phi_block.zk) + 1;

        if (u < w) return 0;

        if (u <= IACBRTX) // can terminate
            return 0;

        assert(u >= w);

        // sieve interval [w,u] fully, then count remaining primes
        std::vector<bool> aux(u - w + 1, 1);

        for (size_t i = 1; ; ++i)
        {
            assert(i < PRIMES.size());
            uint64_t p = PRIMES[i];
            if (p*p > u)
                break;

            // only need to start marking multiples at p^2
            uint64_t j0 = std::max(p*p, p * ceil_div(w, p));
            for (uint64_t j = j0; j <= u; j += p)
            {
                aux[j - w] = 0;
            }
        }

        // add pi(x / pb) where x / pb is in interval
        for (; u >= w; --u)
        {
            // skip loop if u isn't prime
            if (aux[u - w] == 0)
                continue;

            uint64_t y = X / u; // increasing

            if (y >= phi_block.zk) break;

            ++v_defer;
            P2 += phi_block.sum_to(y) + A - 1;
            ++defer;
        }

        #pragma omp atomic
        v += v_defer;

        phi_defer = defer;
        return P2;
    }


    uint64_t primecount(void)
    {
        // Sum terms
        uint64_t S0 = S0_iter();
        uint64_t S1 = 0;
        uint64_t S2 = 0;
        uint64_t P2 = A * (A-1) / 2; // starting sum
        uint64_t v = A;

        std::cout << "S0 = " << S0 << "\n";

        // Init S2 vars
        for (uint64_t b = ASTAR; b < A - 1; ++b)
        {
            uint64_t pb1 = PRIMES[b+1];

            uint64_t tb;
            // hope this is accurate
            if (cbrt(X) <= pb1)     tb = b + 2;
            else if (pb1*pb1 <= Z)  tb = A + 1;
            else                    tb = PRIME_COUNT[X / (pb1*pb1)] + 1;

            S2 += A - (tb - 1); // S2_b
        }

        // block_sum[k][b] = phi(zk - 1, b) - phi(zk1 - 1, b)
        // = sum of indicators of values with first b primes sieved out
        // phi_save[k][b] = phi(zk - 1, b)

        // by indexing k first, consecutive k (for given b) should be far apart
        // to avoid false sharing (in theory)

        // Only take KL-size batches at once to avoid excessive table space
        std::vector<std::vector<uint64_t>> block_sum(KL, std::vector<uint64_t>(A+1, 0));
        auto phi_save = block_sum;

        // deferred counts of phi_save from phi(y,b) calls
        std::vector<std::vector<uint64_t>> S1_defer(KL, std::vector<uint64_t>(ASTAR+1, 0));
        auto S2_defer = block_sum;
        std::vector<uint64_t> P2_defer(KL, 0);

        // Main segmented sieve: blocks Bk = [z_{k-1}, z_k)
        // k batch [k0, k0 + KL)
        // indexing into KL-size table uses k - k0
        for (uint64_t k0 = 1; k0 <= K; k0 += KL)
        {
            std::cout << "Start new K batch at " << k0 << std::endl;

            uint64_t kmax = std::min(k0 + KL, K + 1); // exclusive

            // Dynamic as the block computations are heavily imbalanced
            // for low k
            #pragma omp parallel for schedule(dynamic)
            for (uint64_t k = k0; k < kmax; ++k)
            {
                uint64_t zk1 = zks[k-1];
                uint64_t zk = zks[k];

                // not actually critical, but lock printing makes it nicer
                #pragma omp critical
                {
                    std::cout << "Start block " << k
                         << std::hex << " [0x" << zk1 << ",0x" << zk << ")"
                         << std::dec << std::endl;
                }

                // construct new phi_block with p1, ..., pc already sieved out
                // using phi_yc precomputed (Appendix I)
                // actually not faster than starting at b = 2
                std::vector<bool> ind((zk-zk1)/2);

                for (size_t i = 0; i < ind.size(); ++i)
                    ind[i] = F_C[(zk1 + 2*i) % Q];

                // init new block
                PhiBlock phi_block(ind, zk1, zk);

                // For each b...
                for (uint64_t b = C; b <= A; ++b)
                {
                    uint64_t pb = PRIMES[b];

                    // sieve out p_b for this block
                    if (b > C)
                        phi_block.sieve_out(pb);

                    // update saved block sum for this block
                    block_sum[k-k0][b] = phi_block.sum_to(phi_block.zk - 1);

                    // S1 leaves, b in [C, ASTAR)
                    if (C <= b && b < ASTAR)
                    {
                        #pragma omp atomic
                        S1 += S1_iter(b, phi_block, S1_defer[k-k0][b]);
                    }

                    // S2 leaves, b in [ASTAR, a-1)
                    else if (ASTAR <= b && b < A - 1)
                    {
                        #pragma omp atomic
                        S2 += S2_iter(b, phi_block, S2_defer[k-k0][b]);
                    }

                    // phi2, after sieved out first a primes
                    else if (b == A)
                    {
                        #pragma omp atomic
                        P2 += P2_iter(phi_block, v, P2_defer[k-k0]);
                    }
                }

                std::cout << "End block " << k << std::endl;
            } // end parallel, implicit barrier

            // sum up all deferred phi(y,b) bases sequentially

            for (uint64_t k = k0; k < kmax; ++k)
            {
                for (uint64_t b = C; b <= A; ++b)
                {
                    // accumulate full phi(zk-1,b) from Bk and all previous
                    uint64_t k_prev = (k == k0) ? KL-1 : k-k0-1;
                    uint64_t phi_prev = phi_save[k_prev][b];

                    phi_save[k-k0][b] = phi_prev + block_sum[k-k0][b];

                    if (b < ASTAR)
                        S1    += phi_prev * S1_defer[k-k0][b];
                    else if (b < A-1)
                        S2    += phi_prev * S2_defer[k-k0][b];
                    else if (b == A)
                        P2    += phi_prev * P2_defer[k-k0];
                }
            }
        }

        // Finalize P2
        P2 -= v * (v - 1) /2;

        std::cout << "S1 = " << S1 << std::endl;
        std::cout << "S2 = " << S2 << std::endl;
        std::cout << "P2 = " << P2 << std::endl;

        return S0 + S1 + S2 + A - 1 - P2;
    }
};

