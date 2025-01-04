#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

#include "fenwick_tree.hpp"

using namespace std;

// tuning parameters
#ifdef DEBUG
const int64_t ALPHA = 1;
const int64_t C     = 1;
const int64_t BS    = 50;
#else
const int64_t ALPHA = 6;       // tuning parameter
const int64_t C     = 8;       // precompute phi_c parameter
const int64_t BS    = 1 << 20; // sieve block size
#endif

// global constants
int64_t X;       // Main value to compute pi(X) for
int64_t Q;       // phi_c table size, shouldn't be too large
int64_t Z;       // X^(2/3) / alpha (approx)
int64_t ISQRTX;  // floor(sqrt(X))
int64_t IACBRTX; // floor(alpha cbrt(X))
int64_t a;       // pi(alpha cbrt(X))

// precomputed tables
vector<int64_t> MU_PMIN;     // mu(n) pmin(n) for [1,acbrtx] 
vector<int64_t> PRIMES;      // primes <= acbrtx
vector<int64_t> PRIME_COUNT; // pi(x) over [1,acbrtx] 
vector<int32_t> PHI_C;       // phi(x,c) over [1,Q]

// signum: returns -1, 0, or 1
int sgn(int64_t x)
{
    return (x > 0) - (x < 0);
}

// precompute PRIMES, PRIME_COUNT, and MU_PMIN with a standard sieve to acbrtx
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

// only ints
int64_t cube(int64_t n)
{
    return n * n * n;
}


// contribution of ordinary leaves to phi(x,a)
int64_t S0_compute(void)
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
    cout << "S0 = " << S0 << "\n";
    return S0;
}

// ceil(x/y) for positive integers
int64_t ceil_div(int64_t x, int64_t y)
{
    return x / y + (x % y > 0);
}

// Represents Bk = [zk1, zk) and functions to compute phi(y,b)
struct PhiBlock
{   
    //int64_t k;                // block index
    int64_t bsize;              // block size
    int64_t zk1;                // z_{k-1}, block lower bound (inclusive)
    int64_t zk;                 // z_k, block upper bound (exclusive)
    vector<int64_t> ind;        // 0/1 to track [pmin(y) > pb]
    fenwick_tree phi_sum;       // partial sums
    vector<int64_t> phi_save;   // phi_save(k,b) = phi(zk1-1,b) from prev block
                                // b is only explicitly used here

    // init block at k=1
    // TODO: phi_sum gets overwritten anyway?
    PhiBlock(int64_t bsize, int64_t a) :
        bsize(bsize), 
        zk1(1),
        zk(bsize + 1),
        ind(bsize, 1), 
        phi_sum(ind),
        phi_save(a + 1, 0)
    {}

    void new_block(int64_t k)
    {
        zk1 = bsize * (k-1) + 1;
        zk = bsize * k + 1;
        ind.assign(bsize, 1);
        phi_sum = fenwick_tree(ind);
        
        cout << "block [" << zk1 << "," << zk << ")\n";
    }

    // TODO: save space with odd values only
    int64_t sum_to(int64_t y, int64_t b)
    {
        assert(y >= zk1);
        assert(y < zk);
        return phi_save[b] + phi_sum.sum_to(y - zk1);
    }

    void sieve_out(int64_t pb)
    {
        // sieve out p_b for this block
        for (int64_t j = pb * ceil_div(zk1, pb); j < zk; j += pb)
        {
            if (ind[j-zk1])   // not marked yet
            {
                phi_sum.add_to(j-zk1, -1);
                ind[j-zk1] = 0;
            }
        }        
    }

    void update_save(int64_t b)
    {
        phi_save[b] += phi_sum.sum_to(bsize-1);
    }
    
};

struct Phi2Info
{
    int64_t P2; // phi2 
    int64_t u; // largest p_b not considered yet
    int64_t v; // count number of primes up to sqrt(x)
    int64_t w; // track first integer represented in aux
    vector<bool> aux; // auxiliary sieve to track primes found up to sqrt(x)
    bool done; // flag for not accidentally running when terminated

    // computation of phi2(x,a)
    void alg3(PhiBlock& phi_block)
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
                    {
                        aux[j-w] = false;
                    }
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
        done = true;
        return;
    }
};



int64_t primecount(void)
{
    // Init
    int64_t S1 = 0, S2 = 0;
    int64_t S0 = S0_compute();

    // contribution of special leaves to phi(x,a)

    int64_t astar = 1;
    while(PRIMES[astar+1] * PRIMES[astar+1] <= IACBRTX)
        ++astar;

    cout << "a* = " << astar << endl;

    // init variables used in phi2(x,a) computation
    struct Phi2Info phi2_info = 
    {
        .P2 = a * (a - 1) / 2,
        .u = ISQRTX,
        .v = a,
        .w = phi2_info.u + 1,
        .aux = vector<bool>(),
        .done = false
    };

    // save phi values for summing
    // phi_save(k, b) = phi(z_k - 1, b) from last block
    struct PhiBlock phi_block(BS, a);


    // S2 decreasing d
    vector<int64_t> d(a-1, a);

    // S1 decreasing m
    vector<int64_t> m(astar, IACBRTX);

    //Main segmented sieve: For each interval Bk = [z_{k-1}, z_k)
    for (int64_t k = 1; ; ++k)
    {
        // init new block
        phi_block.new_block(k);
        if (phi_block.zk1 > Z) break;


        // alg1 for each b...
        // start at 0 to sieve out first C primes
        for (int64_t b = 1; b <= a; ++b)
        {
            //cout << "b " << b << endl;
            int64_t pb = PRIMES[b];

            // sieve out p_b for this block
            phi_block.sieve_out(pb);
            
            // S1 leaves
            if (C <= b && b < astar)
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

            // S2, b in [astar, a-1)
            else if (astar <= b && b < a - 1)     
            {
                int64_t pb1 = PRIMES[b+1];
                // d decreasing => y increasing
                
                for (; d[b] > b + 1; --d[b])
                {   
                    int64_t pd = PRIMES[d[b]];
                    int64_t y = X / (pb1 * pd);
                    
                    assert(y >= phi_block.zk1);
                    if (y >= phi_block.zk) break; // handle in future block

                    int64_t phi_yb = 0;
                    
                    //assert(ALPHA > 1); // not needed?

                    // trivial leaves
                    // counting their number tb is not faster here?
                    if (max(X / (pb1*pb1), pb1) < pd && pd <= IACBRTX) 
                    {
                        phi_yb = 1;
                    }
                    // easy leaves
                    else if (max(Z / pb1, pb1) < pd && 
                        pd <= min(X / (pb1*pb1), IACBRTX))
                    {
                        // condition (13) guarantees y <= alpha cbrtx
                        phi_yb = PRIME_COUNT[y] - b + 1;
                    }

                    // hard leaves
                    else 
                    {
                        phi_yb = phi_block.sum_to(y, b);
                    }

                    S2 += phi_yb;
                }
                
            }

            // phi2, sieved out first a primes 
            else if (b == a && !phi2_info.done)
            {
                phi2_info.alg3(phi_block);
            }


            // for next block k
            phi_block.update_save(b);
        }
    }

    // accumulate final results

    int64_t P2 = phi2_info.P2;

    cout << "P2 = " << P2 << endl;
    cout << "S1 = " << S1 << endl;
    cout << "S2 = " << S2 << endl;

    return S0 + S1 + S2 + a - 1 - P2;
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "Usage: ./primecount 1e12\n";
        return 1;
    }

    // setup global constants
    X = atof(argv[1]); // read float from command line
    Z = (cbrt(X) * cbrt(X) / ALPHA);  // approx
    
    // check alpha isn't set too large
    assert(ALPHA <= pow(X, 1/6.));

    ISQRTX = sqrt(X);
    // ensure floating point truncated values are exact floors
    assert(ISQRTX*ISQRTX <= X && (ISQRTX+1)*(ISQRTX+1) > X);

    IACBRTX = ALPHA * cbrt(X);
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

    int64_t pi = primecount();
    cout << pi << endl;
}
