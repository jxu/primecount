#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

#include "fenwick_tree.hpp"

using namespace std;


// signum: returns -1, 0, or 1
int64_t sgn(int64_t x)
{
    return (x > 0) - (x < 0);
}


// ceil(x/y) for positive integers
// TODO: change to round up to multiple
uint64_t ceil_div(uint64_t x, uint64_t y)
{
    return x / y + (x % y > 0);
}


// Represents Bk = [zk1, zk) and functions to compute phi(y,b)
// The physical block size is half the logical size by only storing odd values
// For example, [51, 101) would map to ind [0, 25) via y -> (y-zk1)/2
struct PhiBlock
{   
    size_t bsize;              // (logical) block size
    uint64_t zk1;                // z_{k-1}, block lower bound (inclusive)
    uint64_t zk;                 // z_k, block upper bound (exclusive)
    vector<bool> ind;           // 0/1 to track [pmin(y) > pb]
    fenwick_tree phi_sum;       // data structure for efficient partial sums
    vector<uint64_t> phi_save;   // phi_save(k,b) = phi(zk1-1,b) from prev block
                                // b is only explicitly used here

    // init block at k=1
    // TODO: phi_sum gets overwritten by first new_block anyway?
    PhiBlock(size_t bsize, uint64_t a) :
        bsize(bsize), 
        zk1(1),
        zk(bsize + 1),
        ind(bsize/2, 1), 
        phi_sum(ind),
        phi_save(a + 1, 0)
    {
        assert(bsize % 2 == 0); // requires even size
    }

    void new_block(uint64_t k)
    {
        zk1 = bsize * (k-1) + 1;
        zk = bsize * k + 1;
        ind.assign(bsize/2, 1); 
        phi_sum.reset(ind);
        // does not reset phi_save!
        
        //cout << "block [" << zk1 << "," << zk << ")\n";
    }

    // phi(y,b) compute 
    uint64_t sum_to(uint64_t y, uint64_t b) const
    {
        assert(y >= zk1);
        assert(y < zk);
        return phi_save[b] + phi_sum.sum_to((y - zk1)/2);
    }

    // sieve out p_b for this block
    void sieve_out(uint64_t pb)
    {
        assert(pb > 2);

        uint64_t jstart = pb * ceil_div(zk1, pb);
        if (jstart % 2 == 0) 
            jstart += pb; // ensure odd
        
        for (uint64_t j = jstart; j < zk; j += 2*pb)
        {
            if (ind[(j-zk1)/2])   // not marked yet
            {
                phi_sum.decrease((j-zk1)/2);
                ind[(j-zk1)/2] = 0;
            }
        }        
    }

    void update_save(uint64_t b)
    {
        phi_save[b] += phi_sum.sum_to(bsize/2-1);
    }
};

// Phi2 saved variables
struct Phi2
{
    uint64_t u;    // largest p_b not considered yet
    uint64_t v;    // count number of primes up to sqrt(x)
    uint64_t w;    // track first integer represented in aux
    vector<bool> aux;  // auxiliary sieve to track primes found up to sqrt(x)
    bool done; // flag for not accidentally running when terminated
};

struct Primecount
{
    // tuning parameters
    uint64_t ALPHA;       // tuning parameter
    uint64_t C = 8;       // precompute phi_c parameter
    size_t BS; // sieve block size   

    // global constants
    uint64_t X;       // Main value to compute pi(X) for
    size_t Q;       // phi_c table size, shouldn't be too large
    uint64_t Z;       // X^(2/3) / alpha (approx)
    uint64_t ISQRTX;  // floor(sqrt(X))
    uint64_t IACBRTX; // floor(alpha cbrt(X))
    uint64_t a;       // pi(alpha cbrt(X)) 
    uint64_t astar;   // p_a*^2 = alpha cbrt(X), i.e. a* = pi(sqrt(alpha cbrt X))

    // precomputed tables (assume alpha cbrt X < INT32_MAX)
    vector<int32_t> MU_PMIN;     // mu(n) pmin(n) for [1,acbrtx] 
    vector<uint32_t> PRIMES;      // primes <= acbrtx
    vector<uint32_t> PRIME_COUNT; // pi(x) over [1,acbrtx] 
    vector<uint32_t> PHI_C;       // phi(x,c) over [1,Q]    


    Primecount(uint64_t x, uint64_t alpha, size_t bs) :
        ALPHA(alpha),
        BS(bs),
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
        sieve_mu_prime();
        pre_phi_c();

        a = PRIME_COUNT[IACBRTX];
        cout << "a = " << a << endl;

        assert(PRIMES.size() > (size_t)a + 1); // need p_{a+1}

        // TODO: more efficient?
        astar = 1;
        while(PRIMES[astar+1] * PRIMES[astar+1] <= IACBRTX)
            ++astar;

        cout << "a* = " << astar << endl;
    }

    // precompute PRIMES, PRIME_COUNT, MU_PMIN with a standard sieve to acbrtx
    void sieve_mu_prime(void)
    {
        // Since p_{a+1} may be needed in S2, we introduce fudge factor
        // and hope it's less than the prime gap

        size_t SIEVE_SIZE = IACBRTX + 200;
        
        MU_PMIN.assign(SIEVE_SIZE+1, 1);     // init to 1s
        MU_PMIN[1] = 1000;               // define pmin(1) = +inf
        PRIMES.push_back(1);             // p0 = 1 by convention
        PRIME_COUNT.resize(SIEVE_SIZE+1);    // init values don't matter here

        uint64_t prime_counter = 0;

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
                for (uint64_t i = j*j; i <= SIEVE_SIZE; i += j*j)
                    MU_PMIN[i] = 0;
            }

            PRIME_COUNT[j] = prime_counter;
        }


    }

    // Precompute phi(x,c) 
    void pre_phi_c(void)
    {
        // compute Q as product of first C primes
        Q = 1;
        for (uint64_t i = 1; i <= C; ++i)
            Q *= PRIMES[i];

        PHI_C.resize(Q+1, 1); // index up to Q inclusive
        PHI_C[0] = 0;

        for (uint64_t i = 1; i <= C; ++i)
        {
            uint64_t p = PRIMES[i]; // ith prime, mark multiples as 0
            for (size_t j = p; j <= Q; j += p)
                PHI_C[j] = 0;
        }

        // accumulate
        for (size_t i = 1; i <= Q; ++i)
            PHI_C[i] += PHI_C[i-1];
    }

    // contribution of ordinary leaves to phi(x,a)
    int64_t S0_compute(void)
    {
        int64_t S0 = 0;
        for (uint64_t n = 1; n <= IACBRTX; ++n)
        {
            if ((uint64_t)abs(MU_PMIN[n]) > PRIMES[C])
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
    int64_t S1_iter(uint64_t b, const PhiBlock& phi_block, vector<uint64_t>& m)
    {
        int64_t S1 = 0;
        uint64_t pb1 = PRIMES[b+1];
        // m decreasing
        for (; m[b] * pb1 > IACBRTX; --m[b])
        {
            uint64_t y = X / (m[b] * pb1);

            assert(y >= phi_block.zk1);
            if (y >= phi_block.zk) break;

            if ((uint64_t)abs(MU_PMIN[m[b]]) > pb1)
            {
                S1 -= sgn(MU_PMIN[m[b]]) * phi_block.sum_to(y, b);
            }
        }        
        return S1;
    }

    // Algorithm 2 hell
    void S2_iter(uint64_t b, 
                 const PhiBlock& phi_block, 
                 vector<uint64_t>& d2,
                 vector<char>& t,
                 vector<uint64_t>& S2
    )
    {
        uint64_t pb1 = PRIMES[b+1];
        
        while (d2[b] > b + 1) // step 2, main loop
        {
            uint64_t y = X / (pb1 * PRIMES[d2[b]]);

            if (t[b] == 2) // hard leaves
            {
                if (y >= phi_block.zk) // step 7
                {
                    break; // pause until next block 
                }
                else // step 8 (contribution using phi_block)
                {
                    S2[b] += phi_block.sum_to(y, b);
                    --d2[b];
                    // repeat loop
                }
            } 
            else // t = 0 or 1, easy leaves
            {
                if (y >= IACBRTX) 
                {
                    t[b] = 2;
                    // since t = 2 is set and d2 didn't change, the new loop
                    // will go to step 7
                }                   
                else // step 3/5
                {
                    uint64_t l = PRIME_COUNT[y] - b + 1;
                    
                    if (t[b] == 0) // step 3
                    {
                        // d' + 1 is the smallest d for which (12) holds:
                        // phi(x / (pb1*pd), b) = pi(x / (pb1*pd)) - b + 1
                        uint64_t d_ = PRIME_COUNT[X / (pb1 * PRIMES[b+l])];

                        // step 4
                        if ((PRIMES[d_+1]*PRIMES[d_+1] <= X / pb1) || (d_ <= b))
                        {
                            t[b] = 1; // goto step 6
                        }
                        else // step 5, clustered easy leaves
                        {
                            S2[b] += l * (d2[b] - d_);
                            d2[b] = d_;
                        }
                        
                    }
                    if (t[b] == 1) // t = 1, sparse easy leaves
                    {   
                        // step 6
                        S2[b] += l;
                        --d2[b];
                    }
                }
            }
        }
        return; // terminate
    }

    // Algorithm 3: computation of phi2(x,a)
    void P2_iter(const PhiBlock& phi_block, Phi2& P, uint64_t& P2)
    {
        // step 3 loop (u decrement steps moved here)
        for (; P.u > IACBRTX; --P.u)
        {
            if (P.u < P.w)
            {
                // new aux sieve [w,u] of size IACBRTX+1
                P.w = max((uint64_t)2, P.u - IACBRTX);
                P.aux.assign(P.u - P.w + 1, true);

                for (uint64_t i = 1; ; ++i) 
                {
                    uint64_t p = PRIMES[i];
                    if (p*p > P.u) break;

                    // only need to sieve values starting at p*p within [w,u]
                    uint64_t jstart = max(p*p, p*ceil_div(P.w,p));
                    for (uint64_t j = jstart; j <= P.u; j += p)
                        P.aux[j-P.w] = false;
                }
            }

            // check u to track largest pb not considered yet
            if (P.aux[P.u-P.w]) // prime
            {
                uint64_t y = X / P.u;
                if (y >= phi_block.zk) 
                    return; // finish this block

                // phi(y,a)
                uint64_t phi = phi_block.sum_to(y, a);
                P2 += phi + a - 1;
                ++P.v; // count new prime
            }                  
        }

        // step 3 terminate
        P2 -= P.v * (P.v - 1) / 2;
        P.done = true;
        return;
    }

    uint64_t primecount();
};


uint64_t Primecount::primecount(void)
{
    // Sum accumulators
    int64_t S1 = 0;
    vector<uint64_t> S2(a-1);
    vector<char>  t(a-1);
    uint64_t P2 = a * (a-1) / 2;

    PhiBlock phi_block(BS, a);

    vector<uint64_t> m(astar, IACBRTX); // S1 decreasing m
    vector<uint64_t> d2(a-1); // S2 decreasing d

    Phi2 phi2 = 
    {
        .u = ISQRTX,
        .v = a,
        .w = ISQRTX + 1,
        .aux = vector<bool>(),
        .done = false
    };

    // Ordinary leaves
    int64_t S0 = S0_compute();
    cout << "S0 = " << S0 << "\n";

    // Init S2 vars
    for (uint64_t b = astar; b < a - 1; ++b) 
    {
        uint64_t pb1 = PRIMES[b+1];

        uint64_t tb;
        // hope this is accurate
        if (cbrt(X) <= pb1)     tb = b + 2;
        else if (pb1*pb1 <= Z)  tb = a + 1; 
        else                    tb = PRIME_COUNT[X / (pb1*pb1)] + 1; 

        d2[b] = tb - 1;
        S2[b] = a - d2[b];
        t[b] = 0;
    }

    //Main segmented sieve: For each interval Bk = [z_{k-1}, z_k)
    for (uint64_t k = 1; ; ++k)
    {
        // init new block
        phi_block.new_block(k);
        if (phi_block.zk1 > Z) break;

        // For each b...
        // start at 2 to sieve out odd primes
        assert(C >= 2);
        assert(C <= astar);
        
        for (uint64_t b = 2; b <= a; ++b)
        {
            //cout << "b " << b << endl;
            uint64_t pb = PRIMES[b];

            // sieve out p_b for this block (including <= C)
            phi_block.sieve_out(pb);
          
            // S1 leaves, b in [C, astar)
            if (C <= b && b < astar)
                S1 += S1_iter(b, phi_block, m);

            // S2 leaves, b in [astar, a-1)
            else if (astar <= b && b < a - 1)     
                S2_iter(b, phi_block, d2, t, S2);

            // phi2, after sieved out first a primes 
            else if (b == a && !phi2.done)
                P2_iter(phi_block, phi2, P2);

            // for next block k
            phi_block.update_save(b);
        }
    }

    uint64_t S2_total = 0;
    for (auto x : S2)
        S2_total += x;

    // accumulate final results
    cout << "S1 = " << S1 << endl;
    cout << "S2 = " << S2_total << endl;
    cout << "P2 = " << P2 << endl;

    return S0 + S1 + S2_total + a - 1 - P2;
}

int main(int argc, char* argv[])
{
    if (!(argc == 2 || argc == 4))
    {
        cerr << "Usage: ./primecount X [BLOCKSIZE ALPHA]\n";
        return 1;
    }

    // setup global constants
    // TODO: dynamically adjust constants
    // block size should be O(x^1/3)
    // alpha = O(log^3 x)

    // read float like 1e12 from command line (may not be exact for > 2^53)
    uint64_t X = atof(argv[1]); 
    size_t bs = 1LL << 20; // empirical good block size
    uint64_t alpha = pow(log10(X), 3) / 150; // empirical O(log^3 x) 

    if (argc == 4) // override defaults
    {
        bs = atof(argv[2]);
        alpha = atoi(argv[3]);
    }

    cout << "Computing for X = " << X << endl;
    cout << "Block size = " << bs << endl;
    cout << "Alpha = " << alpha << endl;

    Primecount P(X, alpha, bs);

    cout << P.primecount() << endl;
}
