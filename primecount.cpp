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
const double ALPHA = 40;    // tuning parameter
const int C = 7;            // precompute phi_c parameter
const int64_t BS = 1 << 10; // sieve block size
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
    MU_PMIN.assign(Z+1, 1);     // init to 1s
    PRIMES.push_back(1);        // p0 = 1 by convention
    PRIME_COUNT.resize(Z+1);    // init values don't matter here

    int64_t i, j;
    int64_t prime_counter = 0;

    // sieve of Eratosthenes, modification to keep track of mu sign and pmin
    for (j = 2; j <= Z; ++j) {
        if (MU_PMIN[j] == 1) { // unmarked, so it is prime
            for (i = j; i <= Z; i += j) {
                MU_PMIN[i] = (MU_PMIN[i] == 1) ? -j : -MU_PMIN[i];
            }
        }
    }

    // complete MU_PMIN, compute PRIMES and PRIME_COUNT
    for (j = 2; j <= Z; ++j) {
        if (MU_PMIN[j] == -j) { // prime
            PRIMES.push_back(j); 
            ++prime_counter;

            // mark multiples of p^2 as 0 for mu
            for (i = j*j; i <= Z; i += j*j)
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

    for (int i = 1; i <= C; ++i) {
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

// Algorithm 3: phi2(x,a)
// TODO: use segmented sieving of [1,z]
int64_t phi2_compute(void) 
{
    a = PRIME_COUNT[IACBRTX];
    int64_t v = PRIME_COUNT[ISQRTX];

    int64_t P2 = a*(a-1)/2 - v*(v-1)/2;

    for (int64_t b = a + 1; b <= v; ++b)
        P2 += PRIME_COUNT[X / PRIMES[b]];

    cout << "P2 = " << P2 << "\n";
    return P2;
}

// contribution of ordinary leaves to phi(x,a)
int64_t S0_compute(void)
{
    int64_t S0 = 0;
    for (int64_t n = 1; n <= IACBRTX; ++n) {
        // pmin(1) = +inf
        if (n == 1 || abs(MU_PMIN[n]) > PRIMES[C]) {
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



// Case 2 leaves: Algorithm 2 hell
// TODO: a should be global?
int64_t S2b_compute(const int64_t a, 
                    const int64_t b, 
                    const fenwick_tree& dyn_sieve) 
{
    int64_t pb1 = PRIMES[b+1];
    int64_t xpb12 = X / (pb1*pb1);

    // number of trivial leaves is a + 1 - tb
    int64_t tb;

    if (xpb12 <= pb1) {
        tb = b + 2;
    } else if (xpb12 < IACBRTX) {
        tb = PRIME_COUNT[xpb12] + 1;    
    } else {
        tb = a + 1;
    }

    // step 1
    int64_t d2b = tb - 1; // largest d not considered yet
    int64_t S2b = a - d2b;
    char t = 0;

    while (d2b > b + 1) { // step 2
        int64_t y = X / (pb1 * PRIMES[d2b]);

        if (t == 0) { // step 3, clustered easy leaves
            if (y >= IACBRTX) {
                t = 2;
            } else {
                int64_t l = PRIME_COUNT[y] - b + 1;
                int64_t d_ = PRIME_COUNT[X / (pb1 * PRIMES[b+l])];

                // step 4 
                int64_t pd1 = PRIMES[d_+1];
                if ((pd1 * pd1 <= X / pb1) || d_ <= b) {
                    t = 1; // goto step 6
                } else {
                    S2b += l * (d2b - d_);
                    d2b = d_;
                    continue; // goto 2
                }
            }
        }

        if (t == 1) { // step 5, sparse easy leaves
            if (y >= IACBRTX) {
                t = 2;
            } else {
                int64_t l = PRIME_COUNT[y] - b + 1;
                // step 6
                S2b += l;
                d2b -= 1;
                continue; // goto 2
            }
        } 

        if (t == 2) { // steps 7-9, hard leaves
            S2b += dyn_sieve.sum_to((y-1)/2);
            d2b -= 1;
        }
    }

    return S2b;
}

int64_t primecount(void)
{

    a = PRIME_COUNT[IACBRTX];

    cout << "a = " << a << "\n";

    int64_t P2 = phi2_compute();
    int64_t S0 = S0_compute();

    // contribution of special leaves to phi(x,a)
    int64_t S = 0;
           
    int64_t astar = C;
    while(PRIMES[astar+1] * PRIMES[astar+1] <= IACBRTX) 
        ++astar;

    cout << "a* = " << astar << endl;

    assert(PRIMES[astar]*PRIMES[astar] <= IACBRTX);
    assert(PRIMES[astar+1]*PRIMES[astar+1] > IACBRTX);
    

    // Init alg1 variables for all b
    
    //vector<int64_t> m1(astar, IACBRTX);
    vector<int64_t> S1(astar, 0);

    // save phi values for summing
    // phi_save(b) = phi(z_{k-1}-1, b) from last block up to z_{k-1}
    vector<int64_t> phi_save(astar, 0);

    // For each interval Bk = [z_{k-1}, z_k)
    int64_t zk; 
    for (int64_t k = 1; (zk = BS*k+1) <= Z; ++k) {

        int64_t zk1 = BS*(k-1) + 1;

        cout << "block [" << zk1 << "," << zk << ")\n";

        // init block tree
        vector<int64_t> Bk(BS, 1);
        fenwick_tree phi_block(Bk);

        // sieve out first C primes in block
        for (int64_t i = 1; i <= C; ++i) {
            int64_t p = PRIMES[i];

            for (int64_t j = p * ceil_div(zk1, p); j < zk; j += p) {
                if (Bk[j-zk1]) { // not marked yet
                    phi_block.add_to(j-zk1, -1);
                    Bk[j-zk1] = 0;
                }

            }
        }



        // alg1 for each b...
        for (int64_t b = C; b < astar; ++b) {
            cout << "b " << b << endl;

            // print block
            for (auto x : Bk) cout << x;
            cout << endl;

            for (int i=0; i<BS; ++i) {
                cout << phi_block.sum_to(i) << " ";
            }
            cout << endl;


            cout << "phi_save " << phi_save[b] << endl;
            
            int64_t pb1 = PRIMES[b+1];

            // m decreasing
            for (int64_t m = IACBRTX; m * pb1 > IACBRTX; --m) {

                int64_t y = X / (m * pb1); 

                if (y < zk1 || y >= zk) {
                    continue;
                }
                
                if (abs(MU_PMIN[m]) > pb1) {
                
                    int64_t phi_yb = phi_save[b] + phi_block.sum_to(y-zk1);
                    cout << "y" << y << " phi_block sumto" << phi_block.sum_to(y-zk1) <<  
                    " phi_yb " << phi_yb << endl;

                    S1[b] -= sgn(MU_PMIN[m]) * phi_yb;
                }
                                            
            }

            // for next block k
            phi_save[b] += phi_block.sum_to(BS-1); 

            // sieve out pb1 for this block for next b 

            for (int64_t j = pb1 * ceil_div(zk1, pb1); j < zk; j += pb1) {
                if (Bk[j-zk1]) { // not marked yet
                    phi_block.add_to(j-zk1, -1);
                    Bk[j-zk1] = 0;
                }

            }       
            
        }
  
    }

    // accumulate final results

    int64_t S1_total = 0;
    for (int64_t b = C; b < astar; ++b) {
        S1_total += S1[b];

    }


    S += S1_total;

    cout << "S = " << S << "\n";

    return S0 + S + a - 1 - P2;
}

int main(int argc, char* argv[])
{
    if (argc != 2) {
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
