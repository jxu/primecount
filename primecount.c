// Compile:
// gcc -O3 -Wall -Wextra -o primecount primecount.cpp -lm
// More checks: add -fwrapv

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>

#define BSIZE (1LL << 20) // requires even size
#define PSIZE (BSIZE / 2)

// Globals used in the calculation
// tuning parameters
uint64_t ALPHA;     // tuning parameter
uint64_t C = 8;     // precompute phi_c parameter

// global constants (after setup)
uint64_t X;         // Main value to compute pi(X) for
size_t Q;           // phi_c table size, shouldn't be too large
uint64_t Z;         // X^(2/3) / alpha (approx)
uint64_t ISQRTX;    // floor(sqrt(X))
uint64_t IACBRTX;   // floor(alpha cbrt(X))
uint64_t a;         // pi(alpha cbrt(X)) 
uint64_t astar;     // p_a*^2 = alpha cbrt(X), a* = pi(sqrt(alpha cbrt X))

// precomputed tables (assume alpha cbrt X < INT32_MAX)
int32_t*    MU_PMIN;      // mu(n) pmin(n) for [1,acbrtx] 
uint32_t*   PRIMES;      // primes <= acbrtx
uint32_t    primes_size; 
uint32_t*   PRIME_COUNT; // pi(x) over [1,acbrtx] 
uint32_t*   PHI_C;       // phi(x,c) over [1,Q]    


// Credit: cp-algorithms (Jakob Kobler), e-maxx.ru (Maxim Ivanov)
// customized to save memory 
// only holds unsigned 32-bit values, takes in bools as input
// 1-based tree, indexes [1:len]
typedef uint32_t fenwick_tree[PSIZE+1]; 

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void fenwick_reset(fenwick_tree t, bool* a)
{
    memset(t, 0, PSIZE+1);

    // linear time construction
    for (uint32_t i = 1; i <= PSIZE; ++i)
    {
        t[i] += a[i-1];
        uint32_t r = i + (i & -i);
        if (r <= PSIZE) 
            t[r] += t[i];
    }
}

// 0-based input
uint32_t fenwick_sum_to(const fenwick_tree t, uint32_t r) 
{
    uint32_t s = 0;
    for (++r; r > 0; r -= r & -r)
        s += t[r];
    return s;
}

// 0-based input
void fenwick_decrease(fenwick_tree t, uint32_t i)
{
    for (++i; i <= PSIZE; i += i & -i)
        --t[i];
}

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

typedef struct
{   
    uint64_t         zk1;       // z_{k-1}, block lower bound (inclusive)
    uint64_t         zk;        // z_k, block upper bound (exclusive)
    bool*             ind;       // 0/1 to track [pmin(y) > pb]
    uint32_t*     phi_sum;   // data structure for efficient partial sums
    uint64_t* phi_save;  // phi_save(k,b) = phi(zk1-1,b) from prev block
                                // b is only explicitly used here
} PhiBlock;

PhiBlock* phi_init(uint64_t a)
{
    PhiBlock* p = malloc(sizeof(PhiBlock));
    p->ind = malloc(PSIZE * sizeof(bool));
    p->phi_sum = calloc(PSIZE+1, sizeof(uint32_t));
    p->phi_save = calloc((a+1), sizeof(uint64_t));

    return p;
}

void phi_new_block(PhiBlock* p, uint64_t k)
{
    p->zk1 = BSIZE * (k-1) + 1;
    p->zk = BSIZE * k + 1;

    for (int i = 0; i < PSIZE; ++i)
        p->ind[i] = 1;
    
    fenwick_reset(p->phi_sum, p->ind);
    // does not reset phi_save!
}

// phi(y,b) compute 
uint64_t phi_sum_to(const PhiBlock* p, uint64_t y, uint64_t b)
{
    assert(y >= p->zk1);
    assert(y < p->zk);
    return p->phi_save[b] + fenwick_sum_to(p->phi_sum, (y - p->zk1)/2);
}

// sieve out p_b for this block
void phi_sieve_out(PhiBlock* p, uint64_t pb)
{
    assert(pb > 2);

    uint64_t jstart = pb * ceil_div(p->zk1, pb);
    if (jstart % 2 == 0) 
        jstart += pb; // ensure odd
    
    for (uint64_t j = jstart; j < p->zk; j += 2*pb)
    {
        uint64_t i = (j - p->zk1) / 2;
        if (p->ind[i])   // not marked yet
        {
            fenwick_decrease(p->phi_sum, i);
            p->ind[i] = 0;
        }
    }
}

void phi_update_save(PhiBlock* p, uint64_t b)
{
    p->phi_save[b] += fenwick_sum_to(p->phi_sum, PSIZE-1);
}

// precompute PRIMES, PRIME_COUNT, MU_PMIN with a standard sieve to acbrtx
void sieve_mu_prime(const size_t SIEVE_SIZE)
{
    // init
    MU_PMIN = malloc((SIEVE_SIZE+1) * sizeof(int32_t));
    PRIMES  = malloc((SIEVE_SIZE+1) * sizeof(uint32_t));
    PRIME_COUNT  = malloc((SIEVE_SIZE+1) * sizeof(uint32_t));

    PRIMES[0] = 1;                      // p0 = 1 by convention
    primes_size = 1; 
                    

    for (size_t i = 0; i <= SIEVE_SIZE; ++i)
    {
        MU_PMIN[i] = 1; // init to 1s        
    }

    MU_PMIN[1] = 1000;                  // define pmin(1) = +inf


    // sieve of Eratosthenes, modification to keep track of mu sign and pmin
    for (size_t j = 2; j <= SIEVE_SIZE; ++j)
    {
        if (MU_PMIN[j] == 1)   // unmarked, so it is prime
        {
            for (size_t i = j; i <= SIEVE_SIZE; i += j)
            {
                MU_PMIN[i] = (MU_PMIN[i] == 1) ? -(int32_t)j : -MU_PMIN[i];
            }
        }
    }

    // complete MU_PMIN, compute PRIMES and PRIME_COUNT
    for (size_t j = 2; j <= SIEVE_SIZE; ++j)
    {
        if (MU_PMIN[j] == -(int64_t)j)   // prime
        {
            // push prime
            PRIMES[primes_size++] = j;

            // mark multiples of p^2 as 0 for mu
            for (uint64_t i = j*j; i <= SIEVE_SIZE; i += j*j)
                MU_PMIN[i] = 0;
        }

        PRIME_COUNT[j] = primes_size - 1; // 1-based
    }
}

// Precompute phi(x,c) 
void pre_phi_c(uint64_t C)
{
    // compute Q as product of first C primes
    Q = 1;
    for (uint64_t i = 1; i <= C; ++i)
        Q *= PRIMES[i];

    PHI_C = malloc((Q+1) * sizeof(uint32_t));

    for (uint32_t i = 1; i <= Q; ++i)
        PHI_C[i] = 1;

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

void init(void)
{
    Z = cbrt(X) * cbrt(X) / ALPHA; // approx
    ISQRTX = sqrt(X);
    IACBRTX = ALPHA * cbrt(X);  

    // check alpha isn't set too large
    assert(ALPHA <= pow(X, 1/6.));  

    // hope floating point truncated values are exact floors

    //assert(ISQRTX*ISQRTX <= X);
    //assert((ISQRTX+1)*(ISQRTX+1) > X);
    // may overflow
    //assert(cube(IACBRTX)  <= cube(ALPHA) * X);
    //assert(cube(IACBRTX+1) > cube(ALPHA) * X);

    printf("Z = %ld\n", Z);
    printf("IACBRTX = %ld\n", IACBRTX);
    printf("ISQRTX = %ld\n", ISQRTX);

    // precompute PRIMES, PRIME_COUNT, MU_PMIN

    // Since p_{a+1} may be needed in S2, we introduce fudge factor
    // and hope it's less than the prime gap
    size_t SIEVE_SIZE = IACBRTX + 200;
    sieve_mu_prime(SIEVE_SIZE);


    a = PRIME_COUNT[IACBRTX];

    assert(primes_size > (size_t)a + 1); // need p_{a+1}

    // good enough 
    astar = 1;
    while(PRIMES[astar+1] * PRIMES[astar+1] <= IACBRTX)
        ++astar;

    
    C = MIN(astar, C);
    printf("a = %ld\n", a);
    printf("a* = %ld\n", astar);
    printf("C = %ld\n", C);
    
    assert(C >= 2);
    assert(C <= astar);
    
    pre_phi_c(C);
}


// contribution of ordinary leaves to phi(x,a)
int64_t S0_iter(void)
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
int64_t S1_iter(const PhiBlock* phi_block, 
                const uint64_t b, 
                uint64_t* mb)
{
    int64_t S1 = 0;
    uint64_t pb1 = PRIMES[b+1];
    // m decreasing
    for (; (*mb) * pb1 > IACBRTX; --(*mb))
    {
        uint64_t y = X / ((*mb) * pb1);

        assert(y >= phi_block->zk1);
        if (y >= phi_block->zk) break;

        if ((uint64_t)abs(MU_PMIN[*mb]) > pb1)
        {
            S1 -= sgn(MU_PMIN[*mb]) * phi_sum_to(phi_block, y, b);
        }
    }        
    return S1;
}

// Algorithm 2 hell
int64_t S2_iter(const PhiBlock* phi_block, 
                const uint64_t b, 
                uint64_t* d2b,
                char* t // t[b], different from other tb
)
{
    int64_t S2b = 0;
    uint64_t pb1 = PRIMES[b+1];
    
    while (*d2b > b + 1) // step 2, main loop
    {
        uint64_t y = X / (pb1 * PRIMES[*d2b]);

        if (*t == 2) // hard leaves
        {
            if (y >= phi_block->zk) // step 7
            {
                break; // pause until next block 
            }
            else // step 8 (contribution using phi_block)
            {
                S2b += phi_sum_to(phi_block, y, b);
                --*d2b;
                // repeat loop
            }
        } 
        else // t = 0 or 1, easy leaves
        {
            if (y >= IACBRTX) 
            {
                *t = 2;
                // since t = 2 is set and d2 didn't change, the new loop
                // will go to step 7
            }                   
            else // step 3/5
            {
                uint64_t l = PRIME_COUNT[y] - b + 1;
                
                if (*t == 0) // step 3
                {
                    // d' + 1 is the smallest d for which (12) holds:
                    // phi(x / (pb1*pd), b) = pi(x / (pb1*pd)) - b + 1
                    uint64_t d_ = PRIME_COUNT[X / (pb1 * PRIMES[b+l])];

                    // step 4
                    if ((PRIMES[d_+1]*PRIMES[d_+1] <= X / pb1) || (d_ <= b))
                    {
                        *t = 1; // goto step 6
                    }
                    else // step 5, clustered easy leaves
                    {
                        S2b += l * (*d2b - d_);
                        *d2b = d_;
                    }
                    
                }
                if (*t == 1) // t = 1, sparse easy leaves
                {   
                    // step 6
                    S2b += l;
                    --*d2b;
                }
            }
        }
    }
    return S2b; // terminate
}


// Algorithm 3: computation of phi2(x,a)
int64_t P2_iter(const PhiBlock* phi_block, 
                uint64_t* u,
                uint64_t* v,
                uint64_t* w,
                bool aux[],
                bool* p2done)
{
    int64_t P2 = 0;
    // step 3 loop (u decrement steps moved here)
    for (; *u > IACBRTX; --*u)
    {
        if (*u < *w)
        {
            // new aux sieve [w,u] of size IACBRTX+1
            *w = MAX((uint64_t)2, *u - IACBRTX);

            for (uint64_t i = 0; i < *u - *w + 1; ++i)
                aux[i] = true;
            
            for (uint64_t i = 1; ; ++i) 
            {
                uint64_t p = PRIMES[i];
                if (p*p > *u) break;

                // only need to sieve values starting at p*p within [w,u]
                uint64_t jstart = MAX(p*p, p*ceil_div(*w,p));
                for (uint64_t j = jstart; j <= *u; j += p)
                    aux[j - *w] = false;
            }
        }

        // check u to track largest pb not considered yet
        if (aux[*u - *w]) // prime
        {
            uint64_t y = X / *u;
            if (y >= phi_block->zk) 
                return P2; // finish this block

            // phi(y,a)
            uint64_t phi = phi_sum_to(phi_block, y, a);  
            P2 += phi + a - 1;
            ++*v; // count new prime
        }                  
    }

    // step 3 terminate
    P2 -= *v * (*v - 1) / 2;
    *p2done = true;
    return P2;
}


uint64_t primecount(void)
{
    // Sum accumulators
    int64_t S0 = S0_iter();

    printf("S0 = %ld\n", S0);

    // S1
    int64_t S1 = 0;
    uint64_t* m = malloc(astar * sizeof(uint64_t));

    for (size_t i = 0; i < astar; ++i)
        m[i] = IACBRTX; // S1 decreasing m
    

    // S2 
    int64_t* S2 = calloc(a-1, sizeof(int64_t));
    uint64_t* d2 = calloc(a-1, sizeof(uint64_t)); // S2 decreasing d
    char* t = calloc((a-1), sizeof(char));

    // Phi2
    uint64_t P2 = a * (a-1) / 2; // starting sum
    uint64_t u = ISQRTX;         // largest p_b not considered yet
    uint64_t v = a;              // count number of primes up to sqrt(x)
    uint64_t w = u + 1;          // track first integer represented in aux
    bool* aux = malloc((IACBRTX+1) * sizeof(bool)); // auxiliary sieve to track primes found
    bool p2done = false;         // flag for algorithm terminated

    PhiBlock* phi_block = phi_init(a);


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
        phi_new_block(phi_block, k);

        if (phi_block->zk1 > Z) break;

        // For each b...
        // start at 2 to sieve out odd primes
        
        for (uint64_t b = 2; b <= a; ++b)
        {
            //cout << "b " << b << endl;
            uint64_t pb = PRIMES[b];

            // sieve out p_b for this block (including <= C)
            phi_sieve_out(phi_block, pb);
          
            // S1 leaves, b in [C, astar)
            if (C <= b && b < astar)
                S1 += S1_iter(phi_block, b, &(m[b])); // pass m[b] by ref

            // S2 leaves, b in [astar, a-1)
            else if (astar <= b && b < a - 1)     
                S2[b] += S2_iter(phi_block, b, &(d2[b]), &(t[b]));

            // phi2, after sieved out first a primes 
            else if (b == a && !p2done)
                P2 += P2_iter(phi_block, &u, &v, &w, aux, &p2done);

            // for next block k
            phi_update_save(phi_block, b);
        }
    }

    uint64_t S2_total = 0;
    for (size_t i = 0; i < a-1; ++i)
        S2_total += S2[i];

    // accumulate final results
    printf("S1 = %ld\n", S1);
    printf("S2 = %ld\n", S2_total);
    printf("P2 = %ld\n", P2);
    
    return S0 + S1 + S2_total + a - 1 - P2;
}



int main(int argc, char* argv[])
{
    if (!(argc == 2 || argc == 4))
    {
        fprintf(stderr, "Usage: ./primecount X [BLOCKSIZE ALPHA]\n");
        return 1;
    }

    // read float like 1e12 from command line (may not be exact for > 2^53)
    X = atof(argv[1]); 
    ALPHA = MAX(1., pow(log10(X), 3) / 150); // empirical O(log^3 x) 

    if (argc == 4) // override defaults
    {
        //bs = atof(argv[2]) TODO: remove
        ALPHA = atoi(argv[3]);
    }

    printf("Computing for X = %ld\n", X);
    printf("Block size = %lld\n", BSIZE);
    printf("Alpha = %ld\n", ALPHA);

    init();

    printf("%ld\n", primecount());
}
