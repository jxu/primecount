#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define CEIL_DIV(x,y) (x) / (y) + ((x) % (y) > 0)


typedef uint32_t* fenwick_tree; // light abstraction

const uint32_t MSB_MASK = 1UL << 31;

// tuning parameters
int64_t ALPHA;          // tuning parameter (integer here)
int64_t C;              // precompute phi_c parameter

// global "constants"
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
bool*       F_C;         // f(x,c) = [pmin(x) > p_c], indicators for phi
int64_t*    zks;         // z_k endpoints for interval [1,z]

// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
// Using the MSB for the underlying bool array isn't faster but it's bit fun

// init with input bool array
// result is 1-indexed int array of size len+1, 0-index stores size
fenwick_tree ft_new(bool* ind, uint32_t len)
{
    assert(len + 1 < MSB_MASK);
    // init array with zeros
    fenwick_tree ft = calloc(len + 1, sizeof *ft);

    ft[0] = len;

    // fancy linear time construction
    for (uint32_t i = 1; i <= len; ++i)
    {
        ft[i] = ft[i] + (uint32_t)ind[i-1]; // add input 0/1 to sums
        // should stay with MSB unset
        uint32_t r = i + (i & -i);
        if (r <= len)
            ft[r] += ft[i]; // push forward (ignoring MSB)
        
        ft[i] |= (uint32_t)ind[i-1] << 31; // set MSB with input bool
    }

    return ft;
}

// sum values a[0..r] (0-based)
uint32_t ft_sum_to(fenwick_tree ft, uint32_t r)
{
    assert(r < ft[0]); // TODO: make into function

    uint32_t s = 0;
    // r can go negative here in the weird 0-indexed tree
    for (++r; r > 0; r -= r & -r)
        s += ft[r];
    return s & ~MSB_MASK; // return without MSB
}

// will only decrease if underlying ind[i] = 1
// 0-based input
// ft is passed by value, but underlying t is the same
void ft_try_decrease(fenwick_tree ft, uint32_t i)
{
    assert (i < ft[0]);
    ++i;

    if (ft[i] & MSB_MASK) // if set
    {
        ft[i] &= ~MSB_MASK; // unset
        for (; i <= ft[0]; i += i & -i)
            --ft[i]; // ignore MSB
    }
}

void ft_delete(fenwick_tree ft)
{
    free(ft);
}

// Represents Bk = [zk1, zk) that partition [1, ceil(z)]
// The interval size should be O(iacbrtx) in theory
//
// In the whole computation, Bk processed sequentially for k = 1 to K
// (potentially parallelizable by tracking phi(zk1-1,b) base values used)
// Within each Bk:
//   For b from 1 to a:
//     Sieve out pb
//     Update S1b and S2b

// The physical block size is half the logical size by only storing odd values
// For example, [51, 101) would map to ind [0, 25) via y -> (y-zk1)/2

typedef struct
{
    int64_t          zk1;       // z_{k-1}, block lower bound (inclusive)
    int64_t          zk;        // z_k, block upper bound (exclusive)
    size_t           bsize;     // logical block size
    size_t           psize;     // physical block size

    fenwick_tree     phi_sum;   // data structure for efficient partial sums
} phi_block;


// construct new block without k or b explicitly
// ind should be length psize
phi_block* phi_block_new(bool ind[], size_t zk1, size_t zk)
{
    phi_block* block = malloc(sizeof(phi_block));

    block->zk1 = zk1;
    block->zk = zk;
    block->bsize = zk - zk1;
    block->psize = block->bsize / 2;

    assert(zk1 % 2 == 1);
    assert(block->bsize % 2 == 0);

    block->phi_sum = ft_new(ind, block->psize);

    return block;
}

// translate into actual index into tree
size_t tree_index(const phi_block* B, const int64_t y)
{
    return (y - B->zk1)/2;
}

// sum contribution of this block to phi(y,b)
// = phi(y,b) - phi(zk1 - 1, b) base
int64_t sum_to(const phi_block* B, const int64_t y)
{
    assert(y >= B->zk1);
    assert(y < B->zk);
    return ft_sum_to(B->phi_sum, tree_index(B, y));
}

// sieve out multiples of p_b for this block (including p_b)
void sieve_out(phi_block* block, int64_t pb)
{
    assert(pb > 2); // 2 already sieved by default

    // sieve out pb
    if (block->zk1 <= pb && pb < block->zk)
        ft_try_decrease(block->phi_sum, tree_index(block, pb));

    // now only need to start at pb^2
    // (doesn't really help)
    int64_t j0 = MAX(pb*pb, pb * CEIL_DIV(block->zk1, pb));
    if (j0 % 2 == 0)
        j0 += pb; // ensure odd

    for (int64_t j = j0; j < block->zk; j += 2*pb)
    {
        ft_try_decrease(block->phi_sum, tree_index(block, j));
    }
}

void phi_block_delete(phi_block* block)
{
    ft_delete(block->phi_sum);
    free(block);
}

// signum: returns -1, 0, or 1
int sgn(int64_t x)
{
    return (x > 0) - (x < 0);
}


// convenient pi upper bound when x is beyond iacbrtx table
// (Rosser and Schoenfeld 1962)
double pi_bound(uint64_t x)
{
    if (x <= 1) return 1;
    return 1.25506 * (double)x / log(x);
}

// precompute PRIMES, PRIME_COUNT, MU_PMIN with a standard sieve to acbrtx
void sieve_mu_prime(const size_t SIEVE_SIZE)
{
    // init MU_PMIN to 1s
    MU_PMIN = calloc(SIEVE_SIZE+1, sizeof(int64_t));
    for (size_t i = 0; i < SIEVE_SIZE+1; ++i)
        MU_PMIN[i] = 1;

    MU_PMIN[1] = 1000;                  // define pmin(1) = +inf
    
    int64_t pc = 0; // prime counter (including p0 = 1)
    PRIMES = calloc(SIEVE_SIZE+1, sizeof(int64_t)); // more than necessary
    PRIMES[pc++] = 1; // push p0 = 1 by convention

    PRIME_COUNT = calloc(SIEVE_SIZE+1, sizeof(int64_t));

    // sieve of Eratosthenes, modification to keep track of mu sign and pmin
    for (size_t j = 2; j <= SIEVE_SIZE; ++j)
    {
        if (MU_PMIN[j] == 1)   // unmarked, so it is prime
        {
            for (size_t i = j; i <= SIEVE_SIZE; i += j)
            {
                MU_PMIN[i] = (MU_PMIN[i] == 1) ? -(int64_t)j : -MU_PMIN[i];
            }
        }
    }

    // complete MU_PMIN, compute PRIMES and PRIME_COUNT
    for (size_t j = 2; j <= SIEVE_SIZE; ++j)
    {
        if (MU_PMIN[j] == -(int64_t)j)   // j is prime
        {
            //printf("j prime %ld\n", j);
            PRIMES[pc++] = j; // push j

            // mark multiples of p^2 as 0 for mu
            for (size_t i = j*j; i <= SIEVE_SIZE; i += j*j)
                MU_PMIN[i] = 0;
        }

        PRIME_COUNT[j] = pc - 1; // don't include p0 = 1
    }
}

// Precompute PHI_C = phi(x,c) and F_C = f(n,c)
void pre_phi_c(void)
{
    // compute Q as product of first C primes
    Q = 1;
    for (size_t i = 1; i <= (size_t)C; ++i)
        Q *= PRIMES[i];

    PHI_C = calloc(Q+1, sizeof(int64_t)); // alloc Q+1 zeros
    
    F_C = calloc(Q+1, sizeof *F_C); // Q+1 ones, index up to Q, inclusive
    // F_C[0] = 0
    for (size_t i = 1; i <= (size_t)Q; ++i)
        F_C[i] = 1;

    for (size_t i = 1; i <= (size_t)C; ++i)
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
int64_t phi_yc(const uint64_t y)
{
    return (y / Q) * PHI_C[Q] + PHI_C[y % Q];
}



void primecount_new(uint64_t x, int64_t alpha, int64_t blockmin, int64_t blockmax)
{
    // Init variables
    ALPHA = alpha;
    X = x;
    // Q after sieve
    // hope floating point truncated values are exact floors
    Z = cbrt(x) * cbrt(x) / alpha; // approx
    ISQRTX = sqrt(x);
    IACBRTX = alpha * cbrt(x);
    BLOCKMIN = blockmin;
    BLOCKMAX = blockmax;
   
    C = 8; // default
    KL = 64; // default

    assert(x >= 100); // program not designed for tiny inputs

    // check alpha isn't set too large
    assert(alpha <= pow(x, 1./6));

    printf("Z = %ld\n", Z);
    printf("IACBRTX = %ld\n", IACBRTX);
    printf("ISQRTX = %ld\n", ISQRTX);

    // precompute PRIMES, PRIME_COUNT, MU_PMIN
    // Since p_{a+1} may be needed in S2, leave margin
    size_t SIEVE_SIZE = IACBRTX + 200;
    sieve_mu_prime(SIEVE_SIZE);

    a = PRIME_COUNT[IACBRTX];
    printf("a = %ld\n", a);

    //assert(PRIMES.size() > (size_t)a + 1); // need p_{a+1}

    astar = PRIME_COUNT[(int64_t)(sqrt(ALPHA) * pow(X, 1/6.))];
    printf("a* = %ld\n", astar);

    C = MIN(astar, C);
    printf("C = %ld\n", C);

    // precompute PHI_C tables
    pre_phi_c();

    // create z_k endpoints
    const size_t BSIZE = 1ULL << BLOCKMAX;

    zks = calloc(MAX(2, Z / BLOCKMIN), sizeof(int64_t));
    K = 0;
    zks[K++] = 1; // starting bound
    printf("zks %ld ", zks[0]);

    int64_t zk;
    for (int64_t i = BLOCKMIN; (i < BLOCKMAX); ++i)
    {
        zk = (1ULL << i) + 1;
        zks[K++] = zk;
        printf("%ld ", zk);
        if (zk > Z) break;
    
    }

    if (zk <= Z)
    {
        
        for (size_t i = 1 + (1ULL << BLOCKMAX); i <= Z + BSIZE; i += BSIZE)
        {
            zks[K++] = i;
            printf("%ld ", i);
        }
    }
    --K; // adjust K to have bounds z0 < z1 < ... < zK

    printf("K = %ld\n", K);
}

// contribution of ordinary leaves to phi(x,a)
int64_t S0_iter(void)
{
    int64_t S0 = 0;
    for (size_t n = 1; n <= (size_t)IACBRTX; ++n)
    {
        if (llabs(MU_PMIN[n]) > PRIMES[C])
        {
            S0 += sgn(MU_PMIN[n]) * phi_yc(X / n);
        }
    }

    return S0;
}

// Algorithm 1
int64_t S1_iter(const size_t b, const phi_block* block, int64_t* phi_defer)
{
    int64_t S1 = 0;
    int64_t pb1 = PRIMES[b+1];
    int64_t defer = 0;
    // m decreasing, modified with fixed upper bound
    int64_t mb = MIN((uint64_t)IACBRTX, X / (uint64_t)(block->zk1 * pb1));

    for (; mb * pb1 > IACBRTX; --mb)
    {
        uint64_t y = X / (mb * pb1);

        assert(y >= (uint64_t)(block->zk1));
        if (y >= (uint64_t)(block->zk)) break;

        if (llabs(MU_PMIN[mb]) > pb1)
        {
            S1 -= sgn(MU_PMIN[mb]) * sum_to(block, y);
            defer -= sgn(MU_PMIN[mb]);
        }
    }
    *phi_defer = defer; // write out phi_defer once
    return S1;
}

// Algorithm 2, reworked from leaves formulas
int64_t S2_iter(const int64_t b, const phi_block* block, int64_t* phi_defer)
{
    int64_t S2b = 0;
    uint64_t pb1 = PRIMES[b+1];
    uint64_t zk1 = block->zk1;
    uint64_t zk = block->zk;
    int64_t defer = 0;

    int64_t d = a; // fixed starting point

    // attempt to optimize starting d bound
    // pd <= x / zk1
    d = MIN(d, (int64_t)pi_bound(X / (pb1 * zk1)));
    // non-trivial leaves should satisfy pd <= max(x/pb1^2, pb1)
    d = MIN(d, (int64_t)pi_bound(X / (pb1 * pb1)));

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
        if (y >= (uint64_t)(IACBRTX))
        {
            S2b += sum_to(block, y);
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
    *phi_defer = defer;
    return S2b;
}

// Algorithm 3: computation of phi2(x,a)
// Modification to original algorithm: use aux sieve [u,w] limited to
// y's range in [zk1,zk) rather than acbrtx size
int64_t P2_iter(const phi_block* block, int64_t* v, int64_t* phi_defer)
{
    int64_t P2 = 0;
    int64_t defer = 0;
    int64_t v_defer = 0;

    // maintain aux sieve, u = tracking pb starting at max, y = x / u
    // x/zk < u <= x/zk1
    // iacbrtx < u <= isqrtx
    int64_t u = MIN((uint64_t)ISQRTX, X / block->zk1);
    int64_t w = MAX((uint64_t)IACBRTX, X / block->zk) + 1;

    if (u < w) return 0;

    if (u <= IACBRTX) // can terminate
        return 0;

    assert(u >= w);

    // sieve interval [w,u] fully, then count remaining primes
    bool* aux = calloc(u - w + 1, sizeof(bool));
    for (size_t i = 0; i <= (size_t)(u - w); ++i)
        aux[i] = 1;

    for (size_t i = 1; ; ++i)
    {
        //assert(i < PRIMES.size());
        int64_t p = PRIMES[i];
        if (p*p > u)
            break;

        // only need to start marking multiples at p^2
        int64_t j0 = MAX(p*p, p * CEIL_DIV(w, p));
        for (int64_t j = j0; j <= u; j += p)
            aux[j - w] = 0;
    }

    // add pi(x / pb) where x / pb is in interval
    for (; u >= w; --u)
    {
        // skip loop if u isn't prime
        if (aux[u - w] == 0)
            continue;

        int64_t y = X / u; // increasing

        if (y >= block->zk) break;

        ++v_defer;
        P2 += sum_to(block, y) + a - 1;
        ++defer;
    }

    *v += v_defer;
    *phi_defer = defer;

    free(aux);
    return P2;
}


int64_t primecount(void)
{
    // Sum accumulators
    int64_t S0 = S0_iter();

    printf("S0 = %ld\n", S0);

    // S1
    int64_t S1 = 0;

    // S2
    int64_t* S2 = calloc(a+1, sizeof *S2);
    int64_t S2s = 0;

    // Phi2
    int64_t P2 = a * (a-1) / 2; // starting sum
    int64_t* vs = calloc(K + 1, sizeof *vs);
    int64_t v = a; // counts primes up to sqrt(x) 

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
    // ridiculous 2D indexing macros
    #define IXA(i,j) ((i)*(a+1) + (j))
    #define IXAS(i,j) ((i)*(astar+1) + (j))
    
    // Only take KL-size batches at once to avoid excessive table space
    int64_t* block_sum      = calloc(IXA(KL,0), sizeof(int64_t));
    int64_t* phi_save       = calloc(IXA(KL,0), sizeof(int64_t));
    int64_t* phi_save_prev  = calloc(a+1, sizeof(int64_t));

    // deferred counts of phi_save from phi(y,b) calls
    int64_t* S1_defer       = calloc(IXAS(KL,0), sizeof(int64_t));
    int64_t* S2_defer       = calloc(IXA(KL,0), sizeof(int64_t));
    int64_t* P2_defer       = calloc(KL, sizeof(int64_t));


    // Main segmented sieve: blocks Bk = [z_{k-1}, z_k)
    // k batch [k0, k0 + KL)
    // indexing into KL-size table uses k - k0
    for (int64_t k0 = 1; k0 <= K; k0 += KL)
    {
        printf("Start new K batch at %ld\n", k0);

        int64_t kmax = MIN(k0 + KL, K + 1); // exclusive

        // Dynamic as the block computations are heavily imbalanced
        // for low k
        // Import to schedule dynamically for unbalanced block sizes
        #pragma omp parallel for schedule(dynamic)
        for (int64_t k = k0; k < kmax; ++k)
        {
            int64_t zk1 = zks[k-1];
            int64_t zk = zks[k];
            assert(zk1 < zk);

            // Message may appear broken in multithreading
            printf("Start block %ld [%lx,%lx)\n", k, zk1, zk);

            // construct new phi_block with p1, ..., pc already sieved out
            // using phi_yc precomputed (Appendix I)
            // actually not faster than starting at b = 2
            bool* ind = calloc((zk-zk1)/2, sizeof(bool));

            for (size_t i = 0; i < (size_t)(zk-zk1)/2; ++i)
                ind[i] = F_C[(zk1 + 2*i) % Q];

            // init new block
            phi_block* block = phi_block_new(ind, zk1, zk);

            // For each b...
            for (int64_t b = C; b <= a; ++b)
            {
                int64_t pb = PRIMES[b];

                // sieve out p_b for this block
                if (b > C)
                    sieve_out(block, pb);

                // update saved block sum (k,b) for this block
                block_sum[IXA(k-k0,b)] = sum_to(block, block->zk - 1);

                // S1 leaves, b in [C, astar)
                if ((int64_t)C <= b && b < astar)
                {
                    #pragma omp atomic
                    S1 += S1_iter(b, block, &(S1_defer[IXAS(k-k0,b)]));
                }

                // S2 leaves, b in [astar, a-1)
                else if (astar <= b && b < a - 1)
                {
                    #pragma omp atomic
                    S2[b] += S2_iter(b, block, &(S2_defer[IXA(k-k0,b)]));
                }

                // phi2, after sieved out first a primes
                else if (b == a)
                {
                    #pragma omp atomic
                    P2 += P2_iter(block, &(vs[k]), &(P2_defer[k-k0]));
                }
            }

            // clean up task memory
            phi_block_delete(block);
            free(ind);

            printf("End block %ld\n", k);
        } // end parallel, implicit barrier

        // sum up all deferred phi(y,b) bases sequentially

        for (int64_t k = k0; k < kmax; ++k)
        {
            for (int64_t b = C; b <= a; ++b)
            {
                // accumulate full phi(zk-1,b) from Bk and all previous
                int64_t phi_prev = 
                    (k == k0)
                    ? phi_save_prev[b]
                    : phi_save[IXA(k-k0-1,b)];

                phi_save[IXA(k-k0,b)] = phi_prev + block_sum[IXA(k-k0,b)];

                if (b < astar)
                    S1    += phi_prev * S1_defer[IXAS(k-k0,b)];
                else if (b < a-1)
                    S2[b] += phi_prev * S2_defer[IXA(k-k0,b)];
                else if (b == a)
                    P2    += phi_prev * P2_defer[k-k0];
            }
            v += vs[k];
        }

        // save block_sum for next batch
        for (size_t i = 0; i < (size_t)a+1; ++i)
            phi_save_prev[i] = phi_save[IXA(KL-1,i)];
    }

    // Accumulate final results
    for (int64_t b = 0; b <= a; ++b)
        S2s += S2[b];

    // Finalize P2
    P2 -= v*(v-1)/2;

    printf("v = %ld\n", v);

    free(block_sum);
    free(phi_save);
    free(phi_save_prev);
    free(S1_defer);
    free(S2_defer);
    free(P2_defer);

    printf("S1 = %ld\nS2 = %ld\nP2 = %ld\n", S1, S2s, P2);

    return S0 + S1 + S2s + a - 1 - P2;
}

// Basic tests
// checks ft.sum_to(i) == sum v[0:i]
void check_ft_equal(const fenwick_tree ft, const bool v[])
{
    uint32_t siz = ft[0];
    uint32_t s = 0;
    for (size_t i = 0; i < siz; ++i)
    {
        s += v[i];
        //printf("%d %d\n", ft_sum_to(ft,i), s);
        
        assert(ft_sum_to(ft, i) == s);
    }
}

void test_fenwick_tree()
{
    // example: fenwick tree over array
    bool v1[5] = {1, 1, 0, 1, 1};
    fenwick_tree ft = ft_new(v1, 5);
    check_ft_equal(ft, v1);

    v1[1] = 0;
    ft_try_decrease(ft, 1);
    check_ft_equal(ft, v1);

    v1[1] = 0;
    ft_try_decrease(ft, 1); // should not change
    check_ft_equal(ft, v1);

    v1[4] = 0;
    ft_try_decrease(ft, 4);
    check_ft_equal(ft, v1);

    /*
    // randomized testing
    std::mt19937 rng(1229);
    const int n = 1000;
    std::uniform_int_distribution<> unif1(0, 1);
    std::uniform_int_distribution<> unifn(0, n-1);

    for (int t = 0; t < 10; ++t) // trials
    {
        // randomly fill bool vector
        std::vector<bool> ind(n);
        for (int j = 0; j < n; ++j)
            ind[j] = unif1(rng);

        // init tree from vector
        fenwick_tree ft(ind);

        check_ft_equal(ft, ind);

        // pick random indices to decrease and check FT
        for (int j = 0; j < 100; ++j)
        {
            int x = unifn(rng);
            ind[x] = 0;
            ft.try_decrease(x);
            check_ft_equal(ft, ind);
        }
    }
    */

    printf("Fenwick tree tests passed\n");
}

// Test PhiBlock values without base match a reference
void check_phiyb(phi_block* pb, const int ref[])
{
    for (int64_t i = pb->zk1; i < pb->zk; ++i)
    {
        assert(sum_to(pb, i) == ref[i - pb->zk1]);
    }
}


void test_phi_block()
{
    const int B = 25;
    bool ind[B];

    for (int i = 0; i < B; ++i)
        ind[i] = 1;

    phi_block* pb = phi_block_new(ind, 1, 51);
    // by design, phi block already has b = 1, p_b = 2 sieved out
    // sieved out evens, so remaining are 1,3,5,7,... = 1 mod 2
    const int phi11[50] =
    {
        1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
        11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,
        21,21,22,22,23,23,24,24,25,25
    };

    // sieved out by 3s, so remaining are 1, 5, 7, 11, ... = 1,5 mod 6
    const int phi12[50] =
    {
        1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,7,7,
        7,7,8,8,9,9,9,9,10,10,11,11,11,11,12,12,13,13,13,13,
        14,14,15,15,15,15,16,16,17,17
    };

    check_phiyb(pb, phi11);

    // sieve out b = 2, p_b = 3
    sieve_out(pb, 3);
    check_phiyb(pb, phi12);

    // TODO: automated test
}



int main(int argc, char* argv[])
{
    if (argc == 1) // special test mode
    {
        test_fenwick_tree();
        test_phi_block();
        printf("All tests passed!\n");
        return 0;
    }

    if (!(argc == 2 || argc == 5))
    {
        fprintf(stderr, "Usage: ./primecount X [ALPHA BLOCKMIN BLOCKMAX]\n");
        return 1;
    }

    // setup primecount tuning parameters to pass in

    // read float like 1e12 from command line (may not be exact for > 2^53)
    double Xf = atof(argv[1]);
    if (Xf > (1ll << 53))
    {
        printf("WARNING: atof may not be exact, "
            "and you may need to change parameters for memory\n");
    }

    if (Xf > 1e19)
    {
        fprintf(stderr, "X too big!");
        return 1;
    }

    // convert double to int
    uint64_t X = Xf;
    int64_t alpha = MAX(1., pow(log10(X), 3) / 150); // empirical O(log^3 x)
    int64_t blockmin = 16;
    int64_t blockmax = 24;

    if (argc == 5) // override defaults
    {
        alpha = atoi(argv[2]);
        blockmin = atoi(argv[3]);
        blockmax = atoi(argv[4]);
    }

    printf("Computing for X = %ld\n", X);
    printf("Alpha = %ld\n", alpha);
    printf("BLOCKMIN = %ld\n", blockmin);
    printf("BLOCKMAX = %ld\n", blockmax);

    // main init
    primecount_new(X, alpha, blockmin, blockmax);

    printf("%ld\n", primecount());

    // free globals on exit?
}
