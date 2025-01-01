#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

#include "fenwick_tree.hpp"

using namespace std;

// tuning parameters
#ifdef DEBUG
const double ALPHA = 1;
const int C = 1;
const int64_t BS = 50;
#else
const double ALPHA = 10;    // tuning parameter
const int C = 7;            // precompute phi_c parameter
const int64_t BS = 1 << 24; // sieve block size
#endif

// global constants
int64_t X;       // Main value to compute phi(X) for
int32_t Q;       // phi_c table size, shouldn't be too large
int64_t Z;       // X^(2/3) / alpha (approx)
int64_t ISQRTX;  //
int64_t IACBRTX; // floor(alpha cbrt(X))
int64_t a;       // pi(alpha cbrt(X))

// precomputed tables
vector<int64_t> MU_PMIN;     // mu(n) pmin(n) for [1,Z] TODO
vector<int64_t> PRIMES;      // primes <= Z TODO
vector<int64_t> PRIME_COUNT; // pi(x) over [1,Z] TODO
vector<int32_t> PHI_C;       // phi(x,c) over [1,Q]

// signum: returns -1, 0, or 1
int sgn(int64_t x)
{
    return (x > 0) - (x < 0);
}

// precompute PRIMES, PRIME_COUNT, and MU_PMIN with a standard sieve
// PRIME_COUNT up to sqrt(X) can be computed with a segmented sieve,
// and MU_PMIN only needs values up to IACBRTX, but it's easier/faster to
// store [1,sqrt(X)] directly
void sieve_mu_prime(void)
{
    // for larger values of X, need to subdivide the interval
    MU_PMIN.assign(ISQRTX+1, 1);     // init to 1s
    PRIMES.push_back(1);        // p0 = 1 by convention
    PRIME_COUNT.resize(ISQRTX+1);    // init values don't matter here

    int64_t i, j;
    int64_t prime_counter = 0;

    // sieve of Eratosthenes, modification to keep track of mu sign and pmin
    for (j = 2; j <= ISQRTX; ++j)
    {
        if (MU_PMIN[j] == 1)   // unmarked, so it is prime
        {
            for (i = j; i <= ISQRTX; i += j)
            {
                MU_PMIN[i] = (MU_PMIN[i] == 1) ? -j : -MU_PMIN[i];
            }
        }
    }

    // complete MU_PMIN, compute PRIMES and PRIME_COUNT
    for (j = 2; j <= ISQRTX; ++j)
    {
        if (MU_PMIN[j] == -j)   // prime
        {
            PRIMES.push_back(j);
            ++prime_counter;

            // mark multiples of p^2 as 0 for mu
            for (i = j*j; i <= ISQRTX; i += j*j)
                MU_PMIN[i] = 0;
        }

        PRIME_COUNT[j] = prime_counter;
    }
}

void pre_phi_c(void)
{
    // compute Q as product of first C primes
    Q = 1;
    for (int i = 1; i <= C; ++i)
        Q *= PRIMES[i];

    PHI_C.resize(Q+1, 1); // index up to Q inclusive
    PHI_C[0] = 0;

    for (int i = 1; i <= C; ++i)
    {
        int p = PRIMES[i]; // ith prime, mark multiples as 0
        for (int32_t j = p; j <= Q; j += p)
            PHI_C[j] = 0;
    }

    // accumulate
    for (int32_t i = 1; i <= Q; ++i)
        PHI_C[i] += PHI_C[i-1];
}

// precomputed phi(x,c)
int64_t phi_c(int64_t y)
{
    return (y / Q) * PHI_C[Q] + PHI_C[y % Q];
}

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
        // pmin(1) = +inf
        if (n == 1 || abs(MU_PMIN[n]) > PRIMES[C])
        {
            S0 += sgn(MU_PMIN[n]) * phi_c(X / n);
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

int64_t primecount(void)
{

    a = PRIME_COUNT[IACBRTX];

    cout << "a = " << a << "\n";

    // TODO: compute v on the way instead of sieving to sqrt x
    int64_t v = PRIME_COUNT[ISQRTX];

    // phi2
    int64_t P2 = a*(a-1)/2 - v*(v-1)/2;

    int64_t S0 = S0_compute();

    // contribution of special leaves to phi(x,a)

    int64_t astar = 1;
    while(PRIMES[astar+1] * PRIMES[astar+1] <= IACBRTX)
        ++astar;

    cout << "a* = " << astar << endl;

    assert(PRIMES[astar]*PRIMES[astar] <= IACBRTX);
    assert(PRIMES[astar+1]*PRIMES[astar+1] > IACBRTX);


    // Init
    int64_t S1 = 0, S2 = 0;

    // trivial leaves

    for (int64_t b = astar; b < a - 1; ++b) 
    {
        int64_t tb = 0;
        int64_t pb1 = PRIMES[b+1];

        if (X / (pb1*pb1) <= pb1) 
            tb = b + 2;
        else if (IACBRTX * pb1 * pb1 <= X) 
            tb = a + 1;
        else
            tb = PRIME_COUNT[X / (pb1*pb1)] + 1;


        // broken
        //cout << "b = " << b << " tb = " << tb << endl;
        S2 += a + 1 - tb; 

    }


    // save phi values for summing
    // phi_save(k, b) = phi(z_k - 1, b) from last block
    vector<int64_t> phi_save(a+1, 0);

    // For each interval Bk = [z_{k-1}, z_k)

    for (int64_t k = 1; ; ++k)
    {
        int64_t zk = BS*k + 1;
        int64_t zk1 = BS*(k-1) + 1;
        if (zk1 > Z) break;

        cout << "block [" << zk1 << "," << zk << ")\n";

        // init block tree
        vector<int64_t> Bk(BS, 1);
        fenwick_tree phi_block(Bk);

        // sieve out first C primes in block
        for (int64_t i = 1; i <= C; ++i)
        {
            int64_t p = PRIMES[i];

            for (int64_t j = p * ceil_div(zk1, p); j < zk; j += p)
            {
                if (Bk[j-zk1])   // not marked yet
                {
                    phi_block.add_to(j-zk1, -1);
                    Bk[j-zk1] = 0;
                }
            }
        }


        // alg1 for each b...
        for (int64_t b = C; b <= a; ++b)
        {
            //cout << "b " << b << endl;

            /*
            // print block

            for (auto x : Bk) cout << x;
            cout << endl;

            for (int i=0; i<BS; ++i) {
                cout << phi_block.sum_to(i) << " ";
            }
            cout << endl;
            */

            int64_t pb1 = PRIMES[b+1];



            // S1
            if (b < astar)
            {
                // m decreasing
                for (int64_t m = IACBRTX; m * pb1 > IACBRTX; --m)
                {
                    int64_t y = X / (m * pb1);

                    if (y < zk1)continue;
                    if (y >= zk) break;

                    if (abs(MU_PMIN[m]) > pb1)
                    {
                        int64_t phi_yb = phi_save[b] + phi_block.sum_to(y-zk1);

                        S1 -= sgn(MU_PMIN[m]) * phi_yb;
                    }

                }
            }

            // S2, b in [astar, a-1)
            else if (b < a - 1)     
            {

                // d decreasing
                for (int64_t d = a; d > b + 1; --d)
                {   
                    int64_t pd = PRIMES[d];

                    int64_t y = X / (pb1 * pd);

                    if (y < zk1) continue;
                    if (y >= zk) break;

                    int64_t phi_yb = 0;

                    // shortcuts 
                    //assert(ALPHA > 1);

                    // trivial leaves already handled
                    if (max(X / (pb1*pb1), pb1) < pd && pd <= IACBRTX) 
                    {
                        phi_yb = 0;
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
                        phi_yb = phi_save[b] + phi_block.sum_to(y-zk1);
                    }

                    S2 += phi_yb;
                }
                
            }

            // phi2, sieved out first a primes 
            else if (b == a)
            {
                // renamed b to d here
                for (int64_t d = v; d >= a + 1; --d)
                {
                    int64_t y = X / PRIMES[d];

                    if (y >= zk) break;

                    if (zk1 <= y && y < zk)
                    {
                        int64_t phi = phi_save[a] + phi_block.sum_to(y - zk1);
                        int64_t pi = phi + a - 1;
                        P2 += pi;
                    }
                }

            }

            // for next block k
            phi_save[b] += phi_block.sum_to(BS-1);

            //cout << "phi_save " << phi_save[b] << endl;

            // sieve out p_{b+1} for this block for next b
            for (int64_t j = pb1 * ceil_div(zk1, pb1); j < zk; j += pb1)
            {
                if (Bk[j-zk1])   // not marked yet
                {
                    phi_block.add_to(j-zk1, -1);
                    Bk[j-zk1] = 0;
                }
            }
        }
    }

    // accumulate final results

    cout << "P2 = " << P2 << "\n";
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

    sieve_mu_prime();
    pre_phi_c();

    int64_t pi = primecount();
    cout << pi << '\n';
}
