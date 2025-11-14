#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdint>

using namespace std; // who will stop me?

typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t  i64;

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// Fits exactly in a power of 2 space!
// Using the MSB for the underlying bool array isn't faster but it's bit fun
class fenwick_tree
{
private:
    const u32 MSB_MASK = 1UL << 31;
    vector<u32>   t;   // 0-indexed tree

public:
    // init with input bool array
    fenwick_tree(const vector<bool>& ind) :
        t(ind.size(), 0)
    {
        assert(ind.size() < MSB_MASK);

        // fancy linear time construction
        for (u32 i = 0; i < t.size(); ++i)
        {
            t[i] += ind[i]; // add input 0/1 to sums
            // should stay with MSB unset
            u32 r = i | (i + 1);
            if (r < t.size())
                t[r] += t[i]; // push forward (ignoring MSB)

            t[i] |= ind[i] << 31; // set MSB with input bool
        }
    }

    // sum values a[0..r] (0-based)
    u32 sum_to(u32 r) const
    {
        assert(size_t(r) < t.size());

        u32 s = 0;
        // r can go "negative" here in the weird 0-indexed tree
        for (; r < MSB_MASK; r = (r & (r + 1)) - 1)
            s += t[r];
        return s & ~MSB_MASK; // return without MSB
    }

    // will only decrease if underlying ind[i] = 1
    // 0-based input
    void try_decrease(u32 i)
    {
        assert(i < t.size());

        if (t[i] & MSB_MASK) // if set
        {
            t[i] &= ~MSB_MASK; // unset
            for (; i < t.size(); i |= (i + 1))
                --t[i]; // ignore MSB
        }
    }
};


// ceil(x/y) for positive integers
inline u64 ceil_div(u64 x, u64 y)
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
    u64     zk1;       // z_{k-1}, block lower bound (inclusive)
    u64     zk;        // z_k, block upper bound (exclusive)
    size_t  bsize;     // logical block size
    size_t  psize;     // physical block size

    fenwick_tree    phi_sum;   // data structure for efficient partial sums

    // construct new block without k or b explicitly
    PhiBlock(const vector<bool>& ind, u64 zk1, u64 zk) :
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
    u64 tree_index(const u64 y) const
    {
        return (y - zk1)/2;
    }

    // sum contribution of this block to phi(y,b)
    // = phi(y,b) - phi(zk1 - 1, b) base
    u64 sum_to(u64 y) const
    {
        assert(y >= zk1);
        assert(y < zk);
        return phi_sum.sum_to(tree_index(y));
    }

    // sieve out multiples of p_b for this block (including p_b)
    void sieve_out(u64 pb)
    {
        assert(pb > 2); // 2 already sieved by default

        // sieve out pb
        if (zk1 <= pb && pb < zk)
            phi_sum.try_decrease(tree_index(pb));

        // now only need to start at pb^2
        // (doesn't really help)
        u64 j0 = max(pb*pb, pb * ceil_div(zk1, pb));
        if (j0 % 2 == 0)
            j0 += pb; // ensure odd

        for (u64 j = j0; j < zk; j += 2*pb)
        {
            phi_sum.try_decrease(tree_index(j));
        }
    }
};

// signum: returns -1, 0, or 1
inline int sgn(i64 x)
{
    return (x > 0) - (x < 0);
}


// convenient pi upper bound when x is beyond iacbrtx table
// (Rosser and Schoenfeld 1962)
inline double pi_bound(u64 x)
{
    if (x <= 1) return 1;
    return 1.25506 * x / log(x);
}


class Primecount
{
public:
    // tuning parameters
    u64 ALPHA;    // tuning parameter (integer here)
    u64 C = 8;    // precompute phi_c parameter

    // "global constants"
    u64 X;         // Main value to compute pi(X) for
    u64 Q;         // phi_c table size, shouldn't be too large
    u64 Z;         // X^(2/3) / alpha (approx)
    u64 ISQRTX;    // floor(sqrt(X))
    u64 IACBRTX;   // floor(alpha cbrt(X))
    u64 A;         // pi(alpha cbrt(X))
    u64 ASTAR;     // p_a*^2 = alpha cbrt(X), a* = pi(sqrt(alpha cbrt X))
    u64 BLOCKMIN;  // minimum block size (in bits)
    u64 BLOCKMAX;  // maximum block size (in bits)
    u64 K;         // max k for blocks
    u64 KL = 64;   // how many k to parallelize at once


    // precomputed tables
    vector<i64>        MU_PMIN;     // mu(n) pmin(n) for [2,acbrtx]
    vector<u64>        PRIMES;      // primes <= acbrtx
    vector<u64>        PRIME_COUNT; // pi(x) over [1,acbrtx]
    vector<u64>        PHI_C;       // phi(x,c) over [1,Q]
    vector<bool>       F_C;         // f(x,c) = [pmin(x) > p_c]
    vector<u64>        zks;         // z_k endpoints for interval [1,z]

    Primecount(u64 x, u64 alpha, u64 blockmin, u64 blockmax) :
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
            throw invalid_argument("Program not designed for tiny inputs");

        // check alpha isn't set too large
        if (ALPHA > pow(X, 1/6.))
            throw invalid_argument("Alpha set too large");

        cout << "Z = " << Z << endl;
        cout << "IACBRTX = " << IACBRTX << endl;
        cout << "ISQRTX = " << ISQRTX << endl;

        // precompute PRIMES, PRIME_COUNT, MU_PMIN
        // Since p_{a+1} may be needed in S2, leave margin
        size_t SIEVE_SIZE = IACBRTX + 200;
        sieve_mu_prime(SIEVE_SIZE);

        A = PRIME_COUNT[IACBRTX];
        cout << "a = " << A << endl;

        assert(PRIMES.size() > (size_t)A + 1); // need p_{a+1}

        ASTAR = PRIME_COUNT[u64(sqrt(ALPHA) * pow(X, 1/6.))];
        cout << "a* = " << ASTAR << endl;

        C = min(ASTAR, C);
        cout << "C = " << C << endl;

        // precompute PHI_C tables
        pre_phi_c();

        // create z_k endpoints
        const size_t BSIZE = 1 << BLOCKMAX;
        zks = {1};

        for (u64 i = BLOCKMIN; i < BLOCKMAX; ++i)
            zks.push_back((1 << i) + 1);

        for (u64 i = 1 + BSIZE; i <= Z + BSIZE; i += BSIZE)
            zks.push_back(i);

        K = zks.size() - 1;
        cout << "K = " << K << endl;
    }

    // precompute PRIMES, PRIME_COUNT, MU_PMIN with a standard sieve to acbrtx
    void sieve_mu_prime(size_t SIEVE_SIZE)
    {
        MU_PMIN.assign(SIEVE_SIZE+1, 1);    // init to 1s
        MU_PMIN[1] = 1000;                  // define pmin(1) = +inf
        PRIMES.push_back(1);                // p0 = 1 by convention
        PRIME_COUNT.resize(SIEVE_SIZE+1);   // init values don't matter here

        u64 pc = 0; // prime counter

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
        for (i64 j = 2; j <= i64(SIEVE_SIZE); ++j)
        {
            if (MU_PMIN[j] == -j)   // prime
            {
                PRIMES.push_back(u64(j));
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
            u64 p = PRIMES[i]; // ith prime, mark multiples as 0
            for (u64 j = p; j <= Q; j += p)
                F_C[j] = 0;
        }

        // accumulate
        for (size_t i = 1; i <= Q; ++i)
            PHI_C[i] = F_C[i] + PHI_C[i-1];
    }

    // phi(y,c) can be found quickly from the table
    // y can be 10^19 so needs to be unsigned
    u64 phi_yc(u64 y)
    {
        return (y / Q) * PHI_C[Q] + PHI_C[y % Q];
    }

    // contribution of ordinary leaves to phi(x,a)
    u64 S0_iter(void)
    {
        u64 S0 = 0;
        for (size_t n = 1; n <= size_t(IACBRTX); ++n)
        {
            if (u64(abs(MU_PMIN[n])) > PRIMES[C])
            {
                S0 += sgn(MU_PMIN[n]) * phi_yc(X / n);
            }
        }

        return S0;
    }

    // Algorithm 1
    u64 S1_iter(const u64 b, const PhiBlock& phi_block, u64& phi_defer)
    {
        u64 S1 = 0;
        u64 pb1 = PRIMES[b+1];
        u64 defer = 0;
        // m decreasing, modified with fixed upper bound
        u64 mb = min(IACBRTX, X / (phi_block.zk1 * pb1));

        for (; mb * pb1 > IACBRTX; --mb)
        {
            u64 y = X / (mb * pb1);

            assert(y >= phi_block.zk1);
            if (y >= phi_block.zk) break;

            if (u64(abs(MU_PMIN[mb])) > pb1)
            {
                S1 -= sgn(MU_PMIN[mb]) * phi_block.sum_to(y);
                defer -= sgn(MU_PMIN[mb]);
            }
        }
        phi_defer = defer;
        return S1;
    }

    // Algorithm 2, reworked from leaves formulas
    u64 S2_iter(const u64 b, const PhiBlock& phi_block, u64& phi_defer)
    {
        u64 S2b = 0;
        u64 pb1 = PRIMES[b+1];
        u64 zk1 = phi_block.zk1;
        u64 zk = phi_block.zk;
        u64 defer = 0;

        u64 d = A; // fixed starting point

        // attempt to optimize starting d bound
        // pd <= x / zk1
        d = min(d, u64(pi_bound(X / (pb1 * zk1))));
        // non-trivial leaves should satisfy pd <= max(x/pb1^2, pb1)
        d = min(d, u64(pi_bound(X / (pb1 * pb1))));

        for (; d > b + 1; --d)
        {
            u64 pd = PRIMES[d];
            u64 y = X / (pb1 * pd);

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
    u64 P2_iter(const PhiBlock& phi_block, u64& v, u64& phi_defer)
    {
        u64 P2 = 0;
        u64 defer = 0;
        u64 v_defer = 0;

        // maintain aux sieve, u = tracking pb starting at max, y = x / u
        // x/zk < u <= x/zk1
        // iacbrtx < u <= isqrtx
        u64 u = min(ISQRTX, X / phi_block.zk1);
        u64 w = max(IACBRTX, X / phi_block.zk) + 1;

        if (u < w) return 0;

        if (u <= IACBRTX) // can terminate
            return 0;

        assert(u >= w);

        // sieve interval [w,u] fully, then count remaining primes
        vector<bool> aux(u - w + 1, 1);

        for (size_t i = 1; ; ++i)
        {
            assert(i < PRIMES.size());
            u64 p = PRIMES[i];
            if (p*p > u)
                break;

            // only need to start marking multiples at p^2
            u64 j0 = max(p*p, p * ceil_div(w, p));
            for (u64 j = j0; j <= u; j += p)
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

            u64 y = X / u; // increasing

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


    u64 primecount(void)
    {
        // Sum terms
        u64 S0 = S0_iter();
        u64 S1 = 0;
        u64 S2 = 0;
        u64 P2 = A * (A-1) / 2; // starting sum
        u64 v = A;

        cout << "S0 = " << S0 << "\n";

        // Init S2 vars
        for (u64 b = ASTAR; b < A - 1; ++b)
        {
            u64 pb1 = PRIMES[b+1];

            u64 tb;
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
        vector<vector<u64>> block_sum(KL, vector<u64>(A+1, 0));
        auto phi_save = block_sum;

        // deferred counts of phi_save from phi(y,b) calls
        vector<vector<u64>> S1_defer(KL, vector<u64>(ASTAR+1, 0));
        auto S2_defer = block_sum;
        vector<u64> P2_defer(KL, 0);

        // Main segmented sieve: blocks Bk = [z_{k-1}, z_k)
        // k batch [k0, k0 + KL)
        // indexing into KL-size table uses k - k0
        for (u64 k0 = 1; k0 <= K; k0 += KL)
        {
            cout << "Start new K batch at " << k0 << endl;

            u64 kmax = min(k0 + KL, K + 1); // exclusive

            // Dynamic as the block computations are heavily imbalanced
            // for low k
            #pragma omp parallel for schedule(dynamic)
            for (u64 k = k0; k < kmax; ++k)
            {
                u64 zk1 = zks[k-1];
                u64 zk = zks[k];

                // not actually critical, but lock printing makes it nicer
                #pragma omp critical
                {
                    cout << "Start block " << k
                         << hex << " [0x" << zk1 << ",0x" << zk << ")"
                         << dec << endl;
                }

                // construct new phi_block with p1, ..., pc already sieved out
                // using phi_yc precomputed (Appendix I)
                // actually not faster than starting at b = 2
                vector<bool> ind((zk-zk1)/2);

                for (size_t i = 0; i < ind.size(); ++i)
                    ind[i] = F_C[(zk1 + 2*i) % Q];

                // init new block
                PhiBlock phi_block(ind, zk1, zk);

                // For each b...
                for (u64 b = C; b <= A; ++b)
                {
                    u64 pb = PRIMES[b];

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

                cout << "End block " << k << endl;
            } // end parallel, implicit barrier

            // sum up all deferred phi(y,b) bases sequentially

            for (u64 k = k0; k < kmax; ++k)
            {
                for (u64 b = C; b <= A; ++b)
                {
                    // accumulate full phi(zk-1,b) from Bk and all previous
                    u64 k_prev = (k == k0) ? KL-1 : k-k0-1;
                    u64 phi_prev = phi_save[k_prev][b];

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

        cout << "S1 = " << S1 << endl;
        cout << "S2 = " << S2 << endl;
        cout << "P2 = " << P2 << endl;

        return S0 + S1 + S2 + A - 1 - P2;
    }
};

