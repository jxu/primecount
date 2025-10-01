#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdint>
#include <sstream>

#include "phi_block.hpp"

using namespace std;

// signum: returns -1, 0, or 1
inline int64_t sgn(int64_t x)
{
    return (x > 0) - (x < 0);
}


// convenient pi upper bound when x is beyond iacbrtx table
// (Rosser and Schoenfeld 1962)
inline double pi_bound(double x)
{
    if (x <= 1) return 1;
    return 1.25506 * x / log(x);
}


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
    size_t  Z;         // X^(2/3) / alpha (approx)
    int64_t ISQRTX;    // floor(sqrt(X))
    int64_t IACBRTX;   // floor(alpha cbrt(X))
    int64_t a;         // pi(alpha cbrt(X))
    int64_t astar;     // p_a*^2 = alpha cbrt(X), a* = pi(sqrt(alpha cbrt X))
    size_t  K;         // max k for blocks
 
    // precomputed tables
    vector<int64_t> MU_PMIN;     // mu(n) pmin(n) for [2,acbrtx]
    vector<int64_t> PRIMES;      // primes <= acbrtx
    vector<int64_t> PRIME_COUNT; // pi(x) over [1,acbrtx]
    vector<int64_t> PHI_C;       // phi(x,c) over [1,Q]
    vector<bool>    F_C;         // f(x,c) = [pmin(x) > p_c]
    vector<size_t>  zks;         // z_k endpoints for interval [1,z]

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

        // Since p_{a+1} may be needed in S2, leave margin
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

        // precompute tables
        pre_phi_c(C);


        // create block endpoints for variable size
        const int BLOCK_BITS_MIN = 16;
        const int BLOCK_BITS_MAX = 24;

        // overwrite bsize
        BSIZE = 1 << BLOCK_BITS_MAX;

        zks = {1};

        for (int64_t i = BLOCK_BITS_MIN; i < BLOCK_BITS_MAX; ++i)
        {
            zks.push_back((1 << i) + 1);
        }


        for (size_t i = 1 + BSIZE; i <= Z + BSIZE; i += BSIZE)
        {
            zks.push_back(i);
        }
        
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

        PHI_C.resize(Q+1, 0);
        F_C.resize(Q+1, 1); // index up to Q, inclusive
        F_C[0] = 0;

        for (size_t i = 1; i <= C; ++i)
        {
            int64_t p = PRIMES[i]; // ith prime, mark multiples as 0
            for (size_t j = p; j <= Q; j += p)
                F_C[j] = 0;
        }

        // accumulate
        for (size_t i = 1; i <= Q; ++i)
            PHI_C[i] = F_C[i] + PHI_C[i-1];
    }

    // phi(y,c) can be found quickly from the table
    int64_t phi_yc(int64_t y)
    {
        assert(y >= 0);
        return (y / Q) * PHI_C[Q] + PHI_C[y % Q];
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
                S0 += sgn(MU_PMIN[n]) * phi_yc(y);
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
        // m decreasing, modified with fixed upper bound
        int64_t mb = min(IACBRTX, X / ((int64_t)phi_block.zk1 * pb1));

        for (; mb * pb1 > IACBRTX; --mb)
        {
            size_t y = X / (mb * pb1);

            assert(y >= phi_block.zk1);
            if (y >= phi_block.zk) break;

            if (abs(MU_PMIN[mb]) > pb1)
            {
                S1 -= sgn(MU_PMIN[mb]) * phi_block.sum_to(y);
                phi_defer -= sgn(MU_PMIN[mb]); 
            }
        }
        return S1;
    }

    // Algorithm 2, reworked from leaves formulas
    int64_t S2_iter(
        const int64_t b,
        const PhiBlock& phi_block,
        int64_t& phi_defer)
    {
        int64_t S2b = 0;
        int64_t pb1 = PRIMES[b+1];
        int64_t zk1 = phi_block.zk1;
        int64_t zk = phi_block.zk;

        int64_t d = a; // fixed starting point

        // attempt to optimize starting d bound
        // pd <= x / zk1
        d = min(d, (int64_t)pi_bound(double(X) / (pb1 * zk1)));
        // non-trivial leaves should satisfy pd <= max(x/pb1^2, pb1)
        d = min(d, (int64_t)pi_bound(double(X) / (pb1 * pb1)));
       
        for (; d > b + 1; --d)
        {
            int64_t pd = PRIMES[d];
            // y is increasing as d is decreasing
            int64_t y = X / (pb1 * pd);

            // make sure y is in [zk1,zk)
            if (y < zk1) 
            {
                continue;
            }
            if (y >= zk)
                break;

            // trivial leaves, should be skipped since already counted
            if (max(X / (pb1 * pb1), pb1) < pd)
            {
            }
            // hard leaves
            else if (y >= IACBRTX)
            {
                S2b += phi_block.sum_to(y);
                phi_defer++;
            }

            // easy leaves
            // can't use clustered leaves because of setting d = d' jumps
            else
            {
                // sparse easy leaves
                S2b += PRIME_COUNT[y] - b + 1;
            }
        }

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

        // maintain aux sieve, u = tracking pb starting at max, y = x / u
        // x/zk < u <= x/zk1
        // iacbrtx < u <= isqrtx
        size_t u = min((size_t)ISQRTX, X / phi_block.zk1);
        size_t w = max((size_t)IACBRTX, X / phi_block.zk) + 1;

        if (u < w) return 0;
        
        //cout << "[w,u] " << w << " " << u << endl;

        if (u <= (size_t)IACBRTX) // can terminate
            return 0;

        assert(u >= w);

        // sieve interval [w,u] fully, then count remaining primes
        vector<bool> aux(u - w + 1, 1);

        for (size_t i = 1; ; ++i)
        {
            assert(i < PRIMES.size());
            int64_t p = PRIMES[i];
            if (p*p > u)
                break;
            
            // only need to start marking multiples at p^2
            size_t j0 = max(p*p, p * ceil_div(w, p));
            for (size_t j = j0; j <= u; j += p)
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

            size_t y = X / u; // increasing

            if (y >= phi_block.zk) break;

            ++v;
            P2 += phi_block.sum_to(y) + a - 1;
            phi_defer += 1;
        }

        return P2;
    }


    int64_t primecount(void)
    {
        // Sum accumulators
        int64_t S0 = S0_iter();

        cout << "S0 = " << S0 << "\n";

        // S1
        int64_t S1 = 0;

        // S2
        vector<int64_t> S2(a+1);
        vector<int64_t> d2(a-1); // S2 decreasing d

        // Phi2
        int64_t P2 = a * (a-1) / 2; // starting sum
        vector<int64_t> vs(K + 1, 0);

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
        }



        // block_sum[k][b] = sum of Bk contents after sieving pb
        // = phi(zk - 1, b) - phi(zk1 - 1, b)
        vector<vector<int64_t>> phi_save(K+2, vector<int64_t>(a+1, 0));
        auto block_sum = phi_save;

        // deferred counts of phi_save from phi(y,b) calls
        auto S1_defer = phi_save; 
        auto S2_defer = phi_save;
        auto P2_defer = phi_save;
 
        // Main segmented sieve: blocks Bk = [z_{k-1}, z_k)
        // OpenMP parallelized here! 
        // Dynamic is important because the block computations are 
        // heavily imbalanced.
        #pragma omp parallel for schedule(dynamic)
        for (size_t k = 1; k <= K; ++k)
        {
            // init new block
            size_t zk1 = zks[k-1];
            size_t zk = zks[k];

            // Message may appear broken in multithreading
            ostringstream os;
            os << "Start block " << k << " " 
                << "[" << zk1 << "," << zk << ")" << endl;
            cout << os.str();

            // construct new phi_block with p1, ..., pc already sieved out
            // using phi_yc precomputed (Appendix I)
            // actually not faster than starting at b = 2
            vector<bool> ind((zk-zk1)/2);

            for (size_t i = 0; i < ind.size(); ++i)
            {
                ind[i] = F_C[(zk1 + 2*i) % Q];
            }

            PhiBlock phi_block = PhiBlock(ind, zk1, zk);

            // For each b...
            for (int64_t b = C; b <= a; ++b)
            {
                //cout << "b " << b << endl;
                int64_t pb = PRIMES[b];

                // sieve out p_b for this block
                if ((size_t)b > C)
                    phi_block.sieve_out(pb);

                // update saved block sum for this block
                block_sum[k][b] = phi_block.sum_to(phi_block.zk - 1);

                // S1 leaves, b in [C, astar)
                if ((int64_t)C <= b && b < astar)
                {
                    #pragma omp atomic
                    S1 += S1_iter(b, phi_block, S1_defer[k][b]);
                }

                // S2 leaves, b in [astar, a-1)
                else if (astar <= b && b < a - 1)
                {
                    #pragma omp atomic
                    S2[b] += S2_iter(b, phi_block, S2_defer[k][b]);
                }

                // phi2, after sieved out first a primes
                else if (b == a)
                {
                    #pragma omp atomic
                    P2 += P2_iter(phi_block, vs[k], P2_defer[k][b]); 
                }
            }

            cout << "End block " << k << endl;
        }

        // sum up all deferred phi(y,b) bases sequentially
        int64_t S2s = 0;
        int64_t v = a;

        for (size_t k = 1; k <= K; ++k)
        {
            for (size_t b = C; b <= (size_t)a; ++b)
            {
                // accumulate phi(zk-1,b)
                // phi_save(k,b) = full phi(zk-1,b) from Bk and all previous
                phi_save[k][b] = phi_save[k-1][b] + block_sum[k][b];
                S1    += phi_save[k-1][b] * S1_defer[k][b];
                S2[b] += phi_save[k-1][b] * S2_defer[k][b];
                P2    += phi_save[k-1][b] * P2_defer[k][b];
            }
            v += vs[k];
        }

        for (int64_t b = 0; b <= a; ++b)
        {
            S2s += S2[b];
        }

        // Finalize P2
        P2 -= v*(v-1)/2;

        // accumulate final results
        cout << "S1 = " << S1 << endl;
        cout << "S2 = " << S2s << endl;
        cout << "P2 = " << P2 << endl;

        return S0 + S1 + S2s + a - 1 - P2;
    }
};


