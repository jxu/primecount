#pragma once
#include <math.h>
#include "phi_block.h"


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
    return 1.25506 * (double)x / log(x);
}


typedef struct Primecount
{
    // tuning parameters
    int64_t ALPHA;          // tuning parameter (integer here)
    int64_t C;              // precompute phi_c parameter

    // "global constants"
    uint64_t X;        // Main value to compute pi(X) for
    int64_t Q;         // phi_c table size, shouldn't be too large
    int64_t Z;         // X^(2/3) / alpha (approx)
    int64_t ISQRTX;    // floor(sqrt(X))
    int64_t IACBRTX;   // floor(alpha cbrt(X))
    int64_t a;         // pi(alpha cbrt(X))
    int64_t astar;     // p_a*^2 = alpha cbrt(X), a* = pi(sqrt(alpha cbrt X))
    int64_t BLOCKMIN;  // minimum block size (in bits)
    int64_t BLOCKMAX;  // maximum block size (in bits)
    int64_t K;         // max k for blocks
    int64_t KL;        // how many k to parallelize at once

    // precomputed tables
    int64_t*    MU_PMIN;     // mu(n) pmin(n) for [2,acbrtx]
    int64_t*    PRIMES;      // primes <= acbrtx
    int64_t*    PRIME_COUNT; // pi(x) over [1,acbrtx]
    int64_t*    PHI_C;       // phi(x,c) over [1,Q]
    bool*       F_C;         // f(x,c) = [pmin(x) > p_c]
    int64_t*    zks;         // z_k endpoints for interval [1,z]
} Primecount;
// C = 8
// KL = 64

Primecount* primecount_new(uint64_t x, int64_t alpha, int64_t blockmin, int64_t blockmax)
{
    Primecount* P = malloc(sizeof(P));
    P->ALPHA = alpha;
    P->X = x;
    // Q after sieve
    // hope floating point truncated values are exact floors
    P->Z = (cbrt(x) * cbrt(x) / alpha); // approx
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

    a = PRIME_COUNT[IACBRTX];
    std::cout << "a = " << a << std::endl;

    assert(PRIMES.size() > (size_t)a + 1); // need p_{a+1}

    astar = PRIME_COUNT[int64_t(sqrt(ALPHA) * pow(X, 1/6.))];
    std::cout << "a* = " << astar << std::endl;

    C = std::min(astar, C);
    std::cout << "C = " << C << std::endl;

    // precompute PHI_C tables
    pre_phi_c(C);

    // create z_k endpoints
    const size_t BSIZE = 1 << BLOCKMAX;
    zks = {1};

    for (int64_t i = BLOCKMIN; i < BLOCKMAX; ++i)
        zks.push_back((1 << i) + 1);

    for (size_t i = 1 + BSIZE; i <= Z + BSIZE; i += BSIZE)
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

        int64_t pc = 0; // prime counter

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
        for (size_t j = 2; j <= SIEVE_SIZE; ++j)
        {
            if (MU_PMIN[j] == -(int64_t)j)   // prime
            {
                PRIMES.push_back(j);
                ++pc;

                // mark multiples of p^2 as 0 for mu
                for (size_t i = j*j; i <= SIEVE_SIZE; i += j*j)
                    MU_PMIN[i] = 0;
            }

            PRIME_COUNT[j] = pc;
        }
    }

    // Precompute phi(x,c)
    void pre_phi_c(size_t C)
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
            int64_t p = PRIMES[i]; // ith prime, mark multiples as 0
            for (int64_t j = p; j <= Q; j += p)
                F_C[j] = 0;
        }

        // accumulate
        for (int64_t i = 1; i <= Q; ++i)
            PHI_C[i] = F_C[i] + PHI_C[i-1];
    }

    // phi(y,c) can be found quickly from the table
    int64_t phi_yc(uint64_t y)
    {
        return (y / Q) * PHI_C[Q] + PHI_C[y % Q];
    }

    // contribution of ordinary leaves to phi(x,a)
    int64_t S0_iter(void)
    {
        int64_t S0 = 0;
        for (size_t n = 1; n <= size_t(IACBRTX); ++n)
        {
            if (std::abs(MU_PMIN[n]) > PRIMES[C])
            {
                S0 += sgn(MU_PMIN[n]) * phi_yc(X / n);
            }
        }

        return S0;
    }

    // Algorithm 1
    int64_t S1_iter(
        const size_t b,
        const PhiBlock& phi_block,
        int64_t& phi_defer)
    {
        int64_t S1 = 0;
        int64_t pb1 = PRIMES[b+1];
        int64_t defer = 0;
        // m decreasing, modified with fixed upper bound
        int64_t mb = std::min(uint64_t(IACBRTX), X / uint64_t(phi_block.zk1 * pb1));

        for (; mb * pb1 > IACBRTX; --mb)
        {
            uint64_t y = X / (mb * pb1);

            assert(y >= uint64_t(phi_block.zk1));
            if (y >= uint64_t(phi_block.zk)) break;

            if (std::abs(MU_PMIN[mb]) > pb1)
            {
                S1 -= sgn(MU_PMIN[mb]) * phi_block.sum_to(y);
                defer -= sgn(MU_PMIN[mb]);
            }
        }
        phi_defer = defer;
        return S1;
    }

    // Algorithm 2, reworked from leaves formulas
    int64_t S2_iter(
        const int64_t b,
        const PhiBlock& phi_block,
        int64_t& phi_defer)
    {
        int64_t S2b = 0;
        uint64_t pb1 = PRIMES[b+1];
        uint64_t zk1 = phi_block.zk1;
        uint64_t zk = phi_block.zk;
        int64_t defer = 0;

        int64_t d = a; // fixed starting point

        // attempt to optimize starting d bound
        // pd <= x / zk1
        d = std::min(d, (int64_t)pi_bound(X / (pb1 * zk1)));
        // non-trivial leaves should satisfy pd <= max(x/pb1^2, pb1)
        d = std::min(d, (int64_t)pi_bound(X / (pb1 * pb1)));

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
            if (y >= uint64_t(IACBRTX))
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
    int64_t P2_iter(
        const PhiBlock& phi_block,
        int64_t& v,
        int64_t& phi_defer)
    {
        int64_t P2 = 0;
        int64_t defer = 0;
        int64_t v_defer = 0;

        // maintain aux sieve, u = tracking pb starting at max, y = x / u
        // x/zk < u <= x/zk1
        // iacbrtx < u <= isqrtx
        int64_t u = std::min((uint64_t)ISQRTX, X / phi_block.zk1);
        int64_t w = std::max((uint64_t)IACBRTX, X / phi_block.zk) + 1;

        if (u < w) return 0;

        if (u <= IACBRTX) // can terminate
            return 0;

        assert(u >= w);

        // sieve interval [w,u] fully, then count remaining primes
        vecbool aux(u - w + 1, 1);

        for (size_t i = 1; ; ++i)
        {
            assert(i < PRIMES.size());
            int64_t p = PRIMES[i];
            if (p*p > u)
                break;

            // only need to start marking multiples at p^2
            int64_t j0 = std::max(p*p, p * ceil_div(w, p));
            for (int64_t j = j0; j <= u; j += p)
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

            int64_t y = X / u; // increasing

            if (y >= phi_block.zk) break;

            ++v_defer;
            P2 += phi_block.sum_to(y) + a - 1;
            ++defer;
        }

        v += v_defer;
        phi_defer = defer;
        return P2;
    }


    int64_t primecount(void)
    {
        // Sum accumulators
        int64_t S0 = S0_iter();

        std::cout << "S0 = " << S0 << "\n";

        // S1
        int64_t S1 = 0;

        // S2
        int64_t* S2(a+1);
        int64_t S2s = 0;

        // Phi2
        int64_t P2 = a * (a-1) / 2; // starting sum
        int64_t* vs(K + 1, 0);
        int64_t v = a;

        // Init S2 vars
        for (int64_t b = astar; b < a - 1; ++b)
        {
            int64_t pb1 = PRIMES[b+1];

            int64_t tb;
            // hope this is accurate
            if (cbrt(X) <= pb1)     tb = b + 2;
            else if (pb1*pb1 <= Z)  tb = a + 1;
            else                    tb = PRIME_COUNT[X / (pb1*pb1)] + 1;

            S2[b] = a - (tb - 1);
        }

        // block_sum[k][b] = phi(zk - 1, b) - phi(zk1 - 1, b)
        // = sum of indicators of values with first b primes sieved out
        // phi_save[k][b] = phi(zk - 1, b)

        // by indexing k first, consecutive k (for given b) should be far apart
        // to avoid false sharing (in theory)

        // Only take KL-size batches at once to avoid excessive table space
        std::vector<int64_t*> block_sum(KL, int64_t*(a+1, 0));
        auto phi_save = block_sum;
        int64_t* phi_save_prev(a+1, 0);

        // deferred counts of phi_save from phi(y,b) calls
        std::vector<int64_t*> S1_defer(KL, int64_t*(astar+1, 0));
        auto S2_defer = block_sum;
        int64_t* P2_defer(KL, 0);

        // Main segmented sieve: blocks Bk = [z_{k-1}, z_k)
        // k batch [k0, k0 + KL)
        // indexing into KL-size table uses k - k0
        for (int64_t k0 = 1; k0 <= K; k0 += KL)
        {
            std::cout << "Start new K batch at " << k0 << std::endl;

            int64_t kmax = std::min(k0 + KL, K + 1); // exclusive

            // Dynamic as the block computations are heavily imbalanced
            // for low k
            #pragma omp parallel for schedule(dynamic)
            for (int64_t k = k0; k < kmax; ++k)
            {
                int64_t zk1 = zks[k-1];
                int64_t zk = zks[k];

                // Message may appear broken in multithreading
                std::cout << std::format(
                              "Start block {} [{:#x},{:#x})\n", k, zk1, zk);

                // construct new phi_block with p1, ..., pc already sieved out
                // using phi_yc precomputed (Appendix I)
                // actually not faster than starting at b = 2
                vecbool ind((zk-zk1)/2);

                for (size_t i = 0; i < ind.size(); ++i)
                    ind[i] = F_C[(zk1 + 2*i) % Q];

                // init new block
                PhiBlock phi_block = PhiBlock(ind, zk1, zk);

                // For each b...
                for (int64_t b = C; b <= a; ++b)
                {
                    int64_t pb = PRIMES[b];

                    // sieve out p_b for this block
                    if (b > C)
                        phi_block.sieve_out(pb);

                    // update saved block sum for this block
                    block_sum[k-k0][b] = phi_block.sum_to(phi_block.zk - 1);

                    // S1 leaves, b in [C, astar)
                    if ((int64_t)C <= b && b < astar)
                    {
                        #pragma omp atomic
                        S1 += S1_iter(b, phi_block, S1_defer[k-k0][b]);
                    }

                    // S2 leaves, b in [astar, a-1)
                    else if (astar <= b && b < a - 1)
                    {
                        #pragma omp atomic
                        S2[b] += S2_iter(b, phi_block, S2_defer[k-k0][b]);
                    }

                    // phi2, after sieved out first a primes
                    else if (b == a)
                    {
                        #pragma omp atomic
                        P2 += P2_iter(phi_block, vs[k], P2_defer[k-k0]);
                    }
                }

                std::cout << "End block " << k << std::endl;
            } // end parallel, implicit barrier

            // sum up all deferred phi(y,b) bases sequentially

            for (int64_t k = k0; k < kmax; ++k)
            {
                for (int64_t b = C; b <= a; ++b)
                {
                    // accumulate full phi(zk-1,b) from Bk and all previous
                    int64_t phi_prev = (k == k0)
                                       ? phi_save_prev[b]
                                       : phi_save[k-k0-1][b];

                    phi_save[k-k0][b] = phi_prev + block_sum[k-k0][b];

                    if (b < astar)
                        S1    += phi_prev * S1_defer[k-k0][b];
                    else if (b < a-1)
                        S2[b] += phi_prev * S2_defer[k-k0][b];
                    else if (b == a)
                        P2    += phi_prev * P2_defer[k-k0];
                }
                v += vs[k];
            }

            // save block_sum for next batch
            phi_save_prev = phi_save[KL-1];
        }

        // Accumulate final results
        for (int64_t b = 0; b <= a; ++b)
            S2s += S2[b];

        // Finalize P2
        P2 -= v*(v-1)/2;

        std::cout << "S1 = " << S1 << std::endl;
        std::cout << "S2 = " << S2s << std::endl;
        std::cout << "P2 = " << P2 << std::endl;

        return S0 + S1 + S2s + a - 1 - P2;
    }
};

