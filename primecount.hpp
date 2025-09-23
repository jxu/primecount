// Main implementation as header-only
// Small API for main and testing programs
#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdint>

using namespace std;

const uint32_t MSB_MASK = 1 << 31;

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// only holds 31-bit values, as sign bit it used to track if already marked
// (using the sign-bit saves a little space and locality by not keeping a
// bool array, but not any faster)
class fenwick_tree
{
private:
    size_t           len; // 0-based len
    vector<uint32_t> t;   // 1-based tree, indexes [1:len]
                          // also uses sign bit

public:
    // init array of 1s of length psize
    fenwick_tree(size_t psize) :
        len(psize)
    {
        assert(psize < MSB_MASK);
        reset();
    }

    // reset reconstruct t
    void reset()
    {
        t.assign(len + 1, 0);
        // linear time construction
        for (size_t i = 1; i <= len; ++i)
        {
            t[i] += 1; // start as if big array of 1s
            size_t r = i + (i & -i);
            if (r <= len) 
                t[r] += t[i];
        }
    }

    // sum values a[0..r] (0-based)
    int32_t sum_to(size_t r) const
    {
        assert(r+1 < t.size());
        int32_t s = 0;
        for (++r; r > 0; r -= r & -r)
            s += t[r]; // ignore sign bit
        return s & (~MSB_MASK);
    }

    // will only decrease if t[i] is not already marked
    // 0-based input
    void try_decrease(size_t i)
    {
        assert(i < t.size());
        if (!(t[i] & MSB_MASK)) // not marked
        {
            t[i] |= MSB_MASK; // mark
            for (++i; i <= len; i += i & -i)
                --t[i]; // should work with MSB
        }
    }
};

// signum: returns -1, 0, or 1
inline int64_t sgn(int64_t x)
{
    return (x > 0) - (x < 0);
}


// ceil(x/y) for positive integers
inline int64_t ceil_div(int64_t x, int64_t y)
{
    return x / y + (x % y > 0);
}


// Represents Bk = [zk1, zk) and functions to compute phi(y,b)
// phi(y,b) = count (pos) integers <= y not divisible by the first b primes
// i.e. not p_b-smooth
// The physical block size is half the logical size by only storing odd values
// For example, [51, 101) would map to ind [0, 25) via y -> (y-zk1)/2

class PhiBlock
{
public:
    size_t           bsize;     // logical block size
    size_t           psize;     // physical block size
    size_t           zk1;       // z_{k-1}, block lower bound (inclusive)
    size_t           zk;        // z_k, block upper bound (exclusive)
    fenwick_tree     phi_sum;   // data structure for efficient partial sums

    // phi_save(k,b) = phi(zk1-1,b) from prev block, for implicit b
    // c <= b < a-1
    vector<int64_t> phi_save;

    // init block at k=1
    PhiBlock(size_t a, size_t bsize) :
        bsize(bsize),
        psize(bsize / 2),  // actually only store odd values
        phi_sum(psize),
        phi_save(a + 1, 0)
    {
        assert(bsize % 2 == 0); // only tested even sizes
    }

    void new_block(size_t k)
    {
        // TODO: more flexible zk1 and zk?
        zk1 = bsize * (k-1) + 1;
        zk = bsize * k + 1;
        phi_sum.reset();
        // does not reset phi_save!

        cout << "block [" << zk1 << "," << zk << ")\n";
    }

    // translate into actual index into tree
    size_t tree_index(const size_t y) const
    {
        return (y - zk1)/2;
    }

    // phi(y,b) compute
    int64_t sum_to(size_t y, size_t b) const
    {
        assert(y >= zk1);
        assert(y < zk);
        return phi_save[b] + phi_sum.sum_to(tree_index(y));
    }

    // sieve out p_b for this block
    void sieve_out(size_t pb)
    {
        assert(pb > 2);

        size_t jstart = pb * ceil_div(zk1, pb);
        if (jstart % 2 == 0)
            jstart += pb; // ensure odd

        for (size_t j = jstart; j < zk; j += 2*pb)
        {
            phi_sum.try_decrease(tree_index(j));
        }
    }

    void update_save(size_t b)
    {
        phi_save[b] += phi_sum.sum_to(psize-1);
    }
};

class Primecount
{
public:
    // tuning parameters
    int64_t ALPHA;     // tuning parameter (integer here)
    size_t  C = 8;     // precompute phi_c parameter
    size_t  BSIZE;     // (logical) block size

    // global constants
    int64_t X;         // Main value to compute pi(X) for
    size_t  Q;         // phi_c table size, shouldn't be too large
    int64_t Z;         // X^(2/3) / alpha (approx)
    int64_t ISQRTX;    // floor(sqrt(X))
    int64_t IACBRTX;   // floor(alpha cbrt(X))
    int64_t a;         // pi(alpha cbrt(X))
    int64_t astar;     // p_a*^2 = alpha cbrt(X), a* = pi(sqrt(alpha cbrt X))

    // precomputed tables (assume alpha cbrt X < INT32_MAX)
    vector<int32_t> MU_PMIN;     // mu(n) pmin(n) for [2,acbrtx]
    vector<int32_t> PRIMES;      // primes <= acbrtx
    vector<int32_t> PRIME_COUNT; // pi(x) over [1,acbrtx]
    vector<int32_t> PHI_C;       // phi(x,c) over [1,Q]


    Primecount(int64_t x, int64_t alpha, size_t bsize) :
        ALPHA(alpha),
        BSIZE(bsize),
        X(x),
        // Q after sieve
        Z(cbrt(X) * cbrt(X) / ALPHA), // approx
        ISQRTX(sqrt(X)),
        IACBRTX(ALPHA * cbrt(X))

    {
        // check alpha isn't set too large
        assert(ALPHA <= pow(X, 1/6.));

        // hope floating point truncated values are exact floors

        //assert(ISQRTX*ISQRTX <= X);
        //assert((ISQRTX+1)*(ISQRTX+1) > X);
        // may overflow
        //assert(cube(IACBRTX)  <= cube(ALPHA) * X);
        //assert(cube(IACBRTX+1) > cube(ALPHA) * X);

        cout << "Z = " << Z << endl;
        cout << "IACBRTX = " << IACBRTX << endl;
        cout << "ISQRTX = " << ISQRTX << endl;

        // precompute PRIMES, PRIME_COUNT, MU_PMIN

        // Since p_{a+1} may be needed in S2, we introduce fudge factor
        // and hope it's less than the prime gap
        size_t SIEVE_SIZE = IACBRTX + 200;
        sieve_mu_prime(SIEVE_SIZE);


        a = PRIME_COUNT[IACBRTX];
        cout << "a = " << a << endl;

        assert(PRIMES.size() > (size_t)a + 1); // need p_{a+1}

        // good enough
        astar = 1;
        while(PRIMES[astar+1] * PRIMES[astar+1] <= IACBRTX)
            ++astar;

        cout << "a* = " << astar << endl;
        C = min((size_t)astar, C);
        cout << "C = " << C << endl;

        assert(C >= 2);
        assert((int64_t)C <= astar);

        pre_phi_c(C);
    }

    // precompute PRIMES, PRIME_COUNT, MU_PMIN with a standard sieve to acbrtx
    void sieve_mu_prime(size_t SIEVE_SIZE)
    {
        MU_PMIN.assign(SIEVE_SIZE+1, 1);    // init to 1s
        MU_PMIN[1] = 1000;                  // define pmin(1) = +inf
        PRIMES.push_back(1);                // p0 = 1 by convention
        PRIME_COUNT.resize(SIEVE_SIZE+1);   // init values don't matter here

        int64_t prime_counter = 0;

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
                ++prime_counter;

                // mark multiples of p^2 as 0 for mu
                for (size_t i = j*j; i <= SIEVE_SIZE; i += j*j)
                    MU_PMIN[i] = 0;
            }

            PRIME_COUNT[j] = prime_counter;
        }
    }

    // Precompute phi(x,c)
    void pre_phi_c(size_t C)
    {
        // compute Q as product of first C primes
        Q = 1;
        for (size_t i = 1; i <= C; ++i)
            Q *= PRIMES[i];

        PHI_C.resize(Q+1, 1); // index up to Q inclusive
        PHI_C[0] = 0;

        for (size_t i = 1; i <= C; ++i)
        {
            int64_t p = PRIMES[i]; // ith prime, mark multiples as 0
            for (size_t j = p; j <= Q; j += p)
                PHI_C[j] = 0;
        }

        // accumulate
        for (size_t i = 1; i <= Q; ++i)
            PHI_C[i] += PHI_C[i-1];
    }

    // contribution of ordinary leaves to phi(x,a)
    int64_t S0_iter(void)
    {
        int64_t S0 = 0;
        for (int64_t n = 1; n <= IACBRTX; ++n)
        {
            if (abs(MU_PMIN[n]) > PRIMES[C])
            {
                int64_t y = X / n;
                // use precomputed PHI_C table
                int64_t phi_yc = (y / Q) * PHI_C[Q] + PHI_C[y % Q];
                S0 += sgn(MU_PMIN[n]) * phi_yc;
            }
        }

        return S0;
    }

    // Algorithm 1
    int64_t S1_iter(const size_t b,
                    const PhiBlock& phi_block,
                    int64_t& mb)
    {
        int64_t S1 = 0;
        int64_t pb1 = PRIMES[b+1];
        // m decreasing
        for (; mb * pb1 > IACBRTX; --mb)
        {
            size_t y = X / (mb * pb1);

            assert(y >= phi_block.zk1);
            if (y >= phi_block.zk) break;

            if (abs(MU_PMIN[mb]) > pb1)
            {
                S1 -= sgn(MU_PMIN[mb]) * phi_block.sum_to(y, b);
            }
        }
        return S1;
    }

    // Algorithm 2 hell
    int64_t S2_iter(const int64_t b,
                    const PhiBlock& phi_block,
                    int64_t& d2b,
                    char& t // t[b], different from other tb
                )
    {
        int64_t S2b = 0;
        int64_t pb1 = PRIMES[b+1];

        while (d2b > b + 1) // step 2, main loop
        {
            int64_t y = X / (pb1 * PRIMES[d2b]);

            if (t == 2) // hard leaves
            {
                if ((size_t)y >= phi_block.zk) // step 7
                {
                    break; // pause until next block
                }
                else // step 8 (contribution using phi_block)
                {
                    S2b += phi_block.sum_to(y, b);
                    --d2b;
                    // repeat loop
                }
            }
            else // t = 0 or 1, easy leaves
            {
                if (y >= IACBRTX)
                {
                    t = 2;
                    // since t = 2 is set and d2 didn't change, the new loop
                    // will go to step 7
                }
                else // step 3/5
                {
                    int64_t l = PRIME_COUNT[y] - b + 1;

                    if (t == 0) // step 3
                    {
                        // d' + 1 is the smallest d for which (12) holds:
                        // phi(x / (pb1*pd), b) = pi(x / (pb1*pd)) - b + 1
                        int64_t d_ = PRIME_COUNT[X / (pb1 * PRIMES[b+l])];

                        // step 4
                        if ((PRIMES[d_+1]*PRIMES[d_+1] <= X / pb1) || (d_ <= b))
                        {
                            t = 1; // goto step 6
                        }
                        else // step 5, clustered easy leaves
                        {
                            S2b += l * (d2b - d_);
                            d2b = d_;
                        }

                    }
                    if (t == 1) // t = 1, sparse easy leaves
                    {
                        // step 6
                        S2b += l;
                        --d2b;
                    }
                }
            }
        }
        return S2b; // terminate
    }



    // Algorithm 3: computation of phi2(x,a)
    int64_t P2_iter(const PhiBlock& phi_block,
                    int64_t& u,
                    int64_t& v,
                    int64_t& w,
                    vector<bool>& aux,
                    bool& p2done)
    {
        int64_t P2 = 0;
        // step 3 loop (u decrement steps moved here)
        for (; u > IACBRTX; --u)
        {
            if (u < w)
            {
                // new aux sieve [w,u] of size IACBRTX+1
                w = max(2L, u - IACBRTX);
                aux.assign(u - w + 1, true);

                for (size_t i = 1; ; ++i)
                {
                    int64_t p = PRIMES[i];
                    if (p*p > u) break;

                    // only need to sieve values starting at p*p within [w,u]
                    size_t jstart = max(p*p, p*ceil_div(w,p));
                    for (size_t j = jstart; j <= (size_t)u; j += p)
                        aux[j-w] = false;
                }
            }

            // check u to track largest pb not considered yet
            if (aux[u-w]) // prime
            {
                size_t y = X / u;
                if (y >= phi_block.zk)
                    return P2; // finish this block

                // phi(y,a)
                int64_t phi = phi_block.sum_to(y, a);
                P2 += phi + a - 1;
                ++v; // count new prime
            }
        }

        // step 3 terminate
        P2 -= v * (v - 1) / 2;
        p2done = true;
        return P2;
    }


    int64_t primecount(void)
    {
        // Sum accumulators
        int64_t S0 = S0_iter();

        cout << "S0 = " << S0 << "\n";

        // S1
        int64_t S1 = 0;
        vector<int64_t> m(astar, IACBRTX); // S1 decreasing m

        // S2
        vector<int64_t> S2(a-1);
        vector<int64_t> d2(a-1); // S2 decreasing d
        vector<char>  t(a-1);

        // Phi2
        int64_t P2 = a * (a-1) / 2; // starting sum
        int64_t u = ISQRTX;         // largest p_b not considered yet
        int64_t v = a;              // count number of primes up to sqrt(x)
        int64_t w = u + 1;          // track first integer represented in aux
        vector<bool> aux;            // auxiliary sieve to track primes found
        bool p2done = false;         // flag for algorithm terminated

        PhiBlock phi_block(a, BSIZE);


        // Init S2 vars
        for (int64_t b = astar; b < a - 1; ++b)
        {
            int64_t pb1 = PRIMES[b+1];

            int64_t tb;
            // hope this is accurate
            if (cbrt(X) <= pb1)     tb = b + 2;
            else if (pb1*pb1 <= Z)  tb = a + 1;
            else                    tb = PRIME_COUNT[X / (pb1*pb1)] + 1;

            d2[b] = tb - 1;
            S2[b] = a - d2[b];
            t[b] = 0;
        }

        //Main segmented sieve: For each interval Bk = [z_{k-1}, z_k)
        for (size_t k = 1; ; ++k)
        {
            // init new block
            phi_block.new_block(k);
            if (phi_block.zk1 > (size_t)Z) break;

            // For each b...
            // start at 2 to sieve out odd primes

            for (int64_t b = 2; b <= a; ++b)
            {
                //cout << "b " << b << endl;
                int64_t pb = PRIMES[b];

                // sieve out p_b for this block (including <= C)
                phi_block.sieve_out(pb);

                // S1 leaves, b in [C, astar)
                if ((int64_t)C <= b && b < astar)
                    S1 += S1_iter(b, phi_block, m[b]); // pass m[b] by ref

                // S2 leaves, b in [astar, a-1)
                else if (astar <= b && b < a - 1)
                    S2[b] += S2_iter(b, phi_block, d2[b], t[b]);

                // phi2, after sieved out first a primes
                else if (b == a && !p2done)
                    P2 += P2_iter(phi_block, u, v, w, aux, p2done);

                // for next block k
                phi_block.update_save(b);
            }
        }

        int64_t S2_total = 0;
        for (auto x : S2)
            S2_total += x;

        // accumulate final results
        cout << "S1 = " << S1 << endl;
        cout << "S2 = " << S2_total << endl;
        cout << "P2 = " << P2 << endl;

        return S0 + S1 + S2_total + a - 1 - P2;
    }
};


