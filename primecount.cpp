#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

#include "fenwick_tree.hpp"

using namespace std;


// signum: returns -1, 0, or 1
int sgn(int64_t x)
{
    return (x > 0) - (x < 0);
}


// ceil(x/y) for positive integers
int64_t ceil_div(int64_t x, int64_t y)
{
    return x / y + (x % y > 0);
}

// only ints
int64_t cube(int64_t n)
{
    return n * n * n;
}

// Represents Bk = [zk1, zk) and functions to compute phi(y,b)
// The physical block size is half the logical size by only storing odd values
// For example, [51, 101) would map to ind [0, 25) via y -> (y-zk1)/2
struct PhiBlock
{   
    int64_t bsize;              // (logical) block size
    int64_t zk1;                // z_{k-1}, block lower bound (inclusive)
    int64_t zk;                 // z_k, block upper bound (exclusive)
    vector<bool> ind;           // 0/1 to track [pmin(y) > pb]
    fenwick_tree phi_sum;       // data structure for efficient partial sums
    vector<int64_t> phi_save;   // phi_save(k,b) = phi(zk1-1,b) from prev block
                                // b is only explicitly used here

    // init block at k=1
    // TODO: phi_sum gets overwritten anyway?
    PhiBlock(int64_t bsize, int64_t a) :
        bsize(bsize), 
        zk1(1),
        zk(bsize + 1),
        ind(bsize/2, 1), 
        phi_sum(ind),
        phi_save(a + 1, 0)
    {
        assert(bsize % 2 == 0); // requires even size
    }

    void new_block(int64_t k)
    {
        zk1 = bsize * (k-1) + 1;
        zk = bsize * k + 1;
        ind.assign(bsize/2, 1); 
        phi_sum = fenwick_tree(ind);
        // does not reset phi_save!
        
        cout << "block [" << zk1 << "," << zk << ")\n";
    }

    // phi(y,b) compute 
    int64_t sum_to(int64_t y, int64_t b) const
    {
        assert(y >= zk1);
        assert(y < zk);
        return phi_save[b] + phi_sum.sum_to((y - zk1)/2);
    }

    // sieve out p_b for this block
    void sieve_out(int64_t pb)
    {
        assert(pb > 2);

        int64_t jstart = pb * ceil_div(zk1, pb);
        if (jstart % 2 == 0) 
            jstart += pb; // ensure odd
        
        for (int64_t j = jstart; j < zk; j += 2*pb)
        {
            if (ind[(j-zk1)/2])   // not marked yet
            {
                phi_sum.add_to((j-zk1)/2, -1);
                ind[(j-zk1)/2] = 0;
            }
        }        
    }

    void update_save(int64_t b)
    {
        phi_save[b] += phi_sum.sum_to(bsize/2-1);
    }
};


struct Primecount
{
    // tuning parameters
    int64_t ALPHA;       // tuning parameter
    int64_t C = 8;       // precompute phi_c parameter
    int64_t BS; // sieve block size   

    // global constants
    int64_t X;       // Main value to compute pi(X) for
    int64_t Q;       // phi_c table size, shouldn't be too large
    int64_t Z;       // X^(2/3) / alpha (approx)
    int64_t ISQRTX;  // floor(sqrt(X))
    int64_t IACBRTX; // floor(alpha cbrt(X))
    int64_t a;       // pi(alpha cbrt(X)) 
    int64_t astar;   // p_a*^2 = alpha cbrt(X), i.e. a* = pi(sqrt(alpha cbrt X))

    // precomputed tables
    vector<int64_t> MU_PMIN;     // mu(n) pmin(n) for [1,acbrtx] 
    vector<int64_t> PRIMES;      // primes <= acbrtx
    vector<int64_t> PRIME_COUNT; // pi(x) over [1,acbrtx] 
    vector<int32_t> PHI_C;       // phi(x,c) over [1,Q]    

    // Sum accumulators
    int64_t S0 = 0;
    int64_t S1 = 0;
    vector<int64_t> S2;
    vector<int64_t>  t;
    int64_t P2;

    // Save values between iters
    // save phi values for summing
    // phi_save(k, b) = phi(z_k - 1, b) from last block
    PhiBlock phi_block;

    vector<int64_t> m; // S1 decreasing m
    vector<int64_t> d2; // S2 decreasing d

    // Phi2 saved variables
    int64_t u; // largest p_b not considered yet
    int64_t v; // count number of primes up to sqrt(x)
    int64_t w; // track first integer represented in aux
    vector<bool> aux; // auxiliary sieve to track primes found up to sqrt(x)
    bool phi2done; // flag for not accidentally running when terminated

    Primecount(int64_t x, int64_t alpha, int64_t bs) :
    ALPHA(alpha),
    BS(bs),
    X(x),
    // Q after sieve
    Z(cbrt(X) * cbrt(X) / ALPHA), // approx
    ISQRTX(sqrt(X)),
    IACBRTX(ALPHA * cbrt(X)),

    phi_block(0,0) // empty init

    {
        // check alpha isn't set too large
        assert(ALPHA <= pow(X, 1/6.));

        // ensure floating point truncated values are exact floors
        assert(ISQRTX*ISQRTX <= X && (ISQRTX+1)*(ISQRTX+1) > X);

        assert(cube(IACBRTX)  <= cube(ALPHA) * X &&
               cube(IACBRTX+1) > cube(ALPHA) * X);


        cout << "Computing for X = " << X << endl;
        cout << "Z = " << Z << endl;
        cout << "IACBRTX = " << IACBRTX << endl;
        cout << "ISQRTX = " << ISQRTX << endl;

        
        sieve_mu_prime();
        pre_phi_c();

        a = PRIME_COUNT[IACBRTX];
        cout << "a = " << a << endl;

        astar = 1;
        while(PRIMES[astar+1] * PRIMES[astar+1] <= IACBRTX)
            ++astar;

        cout << "a* = " << astar << endl;

        d2 = vector<int64_t>(a-1, 0);
        m = vector<int64_t>(astar, IACBRTX);

        // init 
        phi_block = PhiBlock(BS, a);

        P2 = a * (a - 1) / 2;
        u = ISQRTX;
        v = a;
        w = u + 1;
        phi2done = false;

        S2.assign(a-1, 0);
        t.assign(a-1, 0);
    }

    // precompute PRIMES, PRIME_COUNT, MU_PMIN with a standard sieve to acbrtx
    void sieve_mu_prime(void)
    {
        // for larger values of X, need to subdivide the interval
        MU_PMIN.assign(IACBRTX+1, 1);     // init to 1s
        MU_PMIN[1] = 1000;               // define pmin(1) = +inf
        PRIMES.push_back(1);             // p0 = 1 by convention
        PRIME_COUNT.resize(IACBRTX+1);    // init values don't matter here

        int64_t i, j;
        int64_t prime_counter = 0;

        // sieve of Eratosthenes, modification to keep track of mu sign and pmin
        for (j = 2; j <= IACBRTX; ++j)
        {
            if (MU_PMIN[j] == 1)   // unmarked, so it is prime
            {
                for (i = j; i <= IACBRTX; i += j)
                {
                    MU_PMIN[i] = (MU_PMIN[i] == 1) ? -j : -MU_PMIN[i];
                }
            }
        }

        // complete MU_PMIN, compute PRIMES and PRIME_COUNT
        for (j = 2; j <= IACBRTX; ++j)
        {
            if (MU_PMIN[j] == -j)   // prime
            {
                PRIMES.push_back(j);
                ++prime_counter;

                // mark multiples of p^2 as 0 for mu
                for (i = j*j; i <= IACBRTX; i += j*j)
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
        for (int64_t i = 1; i <= C; ++i)
            Q *= PRIMES[i];

        PHI_C.resize(Q+1, 1); // index up to Q inclusive
        PHI_C[0] = 0;

        for (int64_t i = 1; i <= C; ++i)
        {
            int64_t p = PRIMES[i]; // ith prime, mark multiples as 0
            for (int64_t j = p; j <= Q; j += p)
                PHI_C[j] = 0;
        }

        // accumulate
        for (int64_t i = 1; i <= Q; ++i)
            PHI_C[i] += PHI_C[i-1];
    }

    // contribution of ordinary leaves to phi(x,a)
    void S0_compute(void)
    {
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
        cout << "S0 = " << S0 << "\n";
    }

    // Algorithm 1
    void S1_iter(int64_t b)
    {
        int64_t pb1 = PRIMES[b+1];
        // m decreasing
        for (; m[b] * pb1 > IACBRTX; --m[b])
        {
            int64_t y = X / (m[b] * pb1);

            assert(y >= phi_block.zk1);
            if (y >= phi_block.zk) break;

            if (abs(MU_PMIN[m[b]]) > pb1)
            {
                S1 -= sgn(MU_PMIN[m[b]]) * phi_block.sum_to(y, b);
            }
        }        
    }

    // Algorithm 2 hell
    void S2_iter(int64_t b, const PhiBlock& phi_block)
    {
        int64_t pb1 = PRIMES[b+1];
        
        while (d2[b] > b + 1) // step 2, main loop
        {
            int64_t y = X / (pb1 * PRIMES[d2[b]]);

            if (t[b] == 2) // hard leaves
            {
                if (y >= phi_block.zk) // step 7
                {
                    return; // pause until next block 
                }
                else // step 8 (contribution using phi_block)
                {
                    S2[b] += phi_block.sum_to(y, b);
                    --d2[b];
                    continue; // repeat loop
                }
            } 
            else // t = 0 or 1, easy leaves
            {
                if (y >= IACBRTX) 
                {
                    t[b] = 2;
                    // since t = 2 is set and d2 didn't change, the new loop
                    // will go to step 7
                    continue;
                }                   
                else // step 3/5
                {
                    int64_t l = PRIME_COUNT[y] - b + 1;
                    
                    if (t[b] == 0) // step 3
                    {
                        // d' + 1 is the smallest d for which (12) holds:
                        // phi(x / (pb1*pd), b) = pi(x / (pb1*pd)) - b + 1
                        int64_t d_ = PRIME_COUNT[X / (pb1 * PRIMES[b+l])];

                        // step 4
                        if ((PRIMES[d_+1]*PRIMES[d_+1] <= X / pb1) || (d_ <= b))
                        {
                            t[b] = 1;
                            // step 6 
                            S2[b] += l;
                            --d2[b];
                            continue;
                        }
                        else // step 5, clustered easy leaves
                        {
                            S2[b] += l * (d2[b] - d_);
                            d2[b] = d_;
                        }
                        
                    }
                    else // t = 1, sparse easy leaves
                    {   
                        // step 6 again
                        S2[b] += l;
                        --d2[b];
                        continue; // not necessary
                    }
                }
            }
        }
        return; // terminate
    }

    // Algorithm 3: computation of phi2(x,a)
    void P2_iter(PhiBlock& phi_block)
    {
        // step 3 loop (u decrement steps moved here)
        for (; u > IACBRTX; --u)
        {
            if (u < w)
            {
                // new aux sieve [w,u] of size IACBRTX+1
                w = max((int64_t)2, u - IACBRTX);
                aux.assign(u - w + 1, true);

                for (int64_t i = 1; ; ++i) 
                {
                    int64_t p = PRIMES[i];
                    if (p*p > u) break;

                    // only need to sieve values starting at p*p within [w,u]
                    for (int64_t j = max(p*p, p*ceil_div(w,p)); j <= u; j += p)
                        aux[j-w] = false;
                }
            }

            // check u to track largest pb not considered yet
            if (aux[u-w]) // prime
            {
                int64_t y = X / u;
                if (y >= phi_block.zk) 
                    return;

                // phi(y,a)
                int64_t phi = phi_block.sum_to(y, a);
                P2 += phi + a - 1;
                ++v; // count new prime
            }                  
        }

        // step 3 terminate
        P2 -= v * (v - 1) / 2;
        phi2done = true;
        return;
    }

    int64_t primecount();
};


int64_t Primecount::primecount(void)
{
    // Ordinary leaves
    S0_compute();

    // Init S2 vars
    for (int64_t b = astar; b < a - 1; ++b) 
    {
        int64_t pb1 = PRIMES[b+1];

        int64_t tb;
        if (X <= cube(pb1))     tb = b + 2;
        else if (pb1*pb1 <= Z)  tb = a + 1; 
        else                    tb = PRIME_COUNT[X / (pb1*pb1)] + 1; 

        d2[b] = tb - 1;
        S2[b] = a - d2[b];
        t[b] = 0;
    }

    //Main segmented sieve: For each interval Bk = [z_{k-1}, z_k)
    for (int64_t k = 1; ; ++k)
    {
        // init new block
        phi_block.new_block(k);
        if (phi_block.zk1 > Z) break;

        // For each b...
        // start at 2 to sieve out odd primes
        assert(C >= 2);
        assert(C < astar);
        
        for (int64_t b = 2; b <= a; ++b)
        {
            //cout << "b " << b << endl;
            int64_t pb = PRIMES[b];

            // sieve out p_b for this block (including <= C)
            phi_block.sieve_out(pb);
          
            // S1 leaves, b in [C, astar)
            if (C <= b && b < astar)
                S1_iter(b);

            // S2 leaves, b in [astar, a-1)
            else if (astar <= b && b < a - 1)     
                S2_iter(b, phi_block);

            // phi2, after sieved out first a primes 
            else if (b == a && !phi2done)
                P2_iter(phi_block);

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

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "Usage: ./primecount 1e12\n";
        return 1;
    }

    // setup global constants
    int64_t X = atof(argv[1]); // read float from command line
    int64_t bs = 1 << 16;
    int64_t alpha = 6;

    Primecount P(X, alpha, bs);

    cout << P.primecount() << endl;
}
