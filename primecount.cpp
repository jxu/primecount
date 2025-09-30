#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdint>

using namespace std;


// Credit: cp-algorithms (Jakob Kogler), e-maxx.ru (Maxim Ivanov)
// customized to save memory by only operating over a bit array (0/1 input)
class fenwick_tree
{
private:
    size_t           len; // 0-based len
    vector<bool>     ind; // underlying array
    vector<uint32_t> t;   // 1-based tree, indexes [1:len]

public:
    // init with input bool array
    fenwick_tree(const vector<bool>& ind) :
        len(ind.size()),
        ind(ind),
        t(len + 1, 0)
    {
        // linear time construction
        for (size_t i = 1; i <= len; ++i)
        {
            t[i] += ind[i-1];
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
            s += t[r];
        return s;
    }

    // will only decrease if ind[i] = 1
    // 0-based input
    void try_decrease(size_t i)
    {
        assert(i < t.size());
        if (ind[i])
        {
            ind[i] = 0;
            for (++i; i <= len; i += i & -i)
                --t[i];
        }
    }
};

// signum: returns -1, 0, or 1
int64_t sgn(int64_t x)
{
    return (x > 0) - (x < 0);
}


// ceil(x/y) for positive integers
int64_t ceil_div(int64_t x, int64_t y)
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
    size_t           zk1;       // z_{k-1}, block lower bound (inclusive)
    size_t           zk;        // z_k, block upper bound (exclusive)
    size_t           bsize;     // logical block size
    size_t           psize;     // physical block size
 
    fenwick_tree     phi_sum;   // data structure for efficient partial sums

    // construct new block without k or b explicitly
    PhiBlock(const vector<bool>& ind, size_t zk1, size_t zk) :
        zk1(zk1),
        zk(zk),
        bsize(zk - zk1),
        psize(bsize / 2),
        phi_sum(ind)
    {
        assert(ind.size() == psize);
        assert(bsize % 2 == 0);
        cout << "block [" << zk1 << "," << zk << ")\n";
    }

    // translate into actual index into tree
    size_t tree_index(const size_t y) const
    {
        return (y - zk1)/2;
    }

    // sum contribution of this block to phi(y,b)
    // = phi(y,b) - phi(zk1 - 1, b) base
    int64_t sum_to(size_t y) const
    {
        assert(y >= zk1);
        assert(y < zk);
        return phi_sum.sum_to(tree_index(y));
    }

    // sieve out multiples of p_b for this block (including p_b)
    void sieve_out(size_t pb)
    {
        assert(pb > 2); // 2 already sieved by default

        size_t jstart = pb * ceil_div(zk1, pb);
        if (jstart % 2 == 0)
            jstart += pb; // ensure odd

        for (size_t j = jstart; j < zk; j += 2*pb)
        {
            phi_sum.try_decrease(tree_index(j));
        }
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
    size_t  Z;         // X^(2/3) / alpha (approx)
    int64_t ISQRTX;    // floor(sqrt(X))
    int64_t IACBRTX;   // floor(alpha cbrt(X))
    int64_t a;         // pi(alpha cbrt(X))
    int64_t astar;     // p_a*^2 = alpha cbrt(X), a* = pi(sqrt(alpha cbrt X))

    // precomputed tables
    vector<int64_t> MU_PMIN;     // mu(n) pmin(n) for [2,acbrtx]
    vector<int64_t> PRIMES;      // primes <= acbrtx
    vector<int64_t> PRIME_COUNT; // pi(x) over [1,acbrtx]
    vector<int64_t> PHI_C;       // phi(x,c) over [1,Q]
    vector<bool>    F_C;         // f(x,c) = [pmin(x) > p_c]


    size_t K; // max k for blocks


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

    // phi(y,b) = sum over n <= y of [p_min(n) > p_b]
    // i.e. count not p_b-smooth
    //int64_t phiyb(const PhiBlock& phi_block, int64_t y, int64_t b) const
    //{
     //   return phi_block.sum_to(y) + phi_save[b];
    //}

    // Algorithm 1
    int64_t S1_iter(
        const size_t b,
        const PhiBlock& phi_block,
        int64_t mb,
        int64_t& phi_defer)
    {
        int64_t S1 = 0;
        int64_t pb1 = PRIMES[b+1];
        // m decreasing, modified with fixed upper bound
        mb = min(mb, X / (int64_t)(phi_block.zk1 * pb1));

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

    // Algorithm 2 hell
    int64_t S2_iter(
        const int64_t b,
        const PhiBlock& phi_block,
        int64_t d2b,
        char t, // t[b], different from other tb
        int64_t& phi_defer)
    {
        int64_t S2b = 0;
        assert(b+1 < PRIMES.size());
        int64_t pb1 = PRIMES[b+1];
        size_t zk1 = phi_block.zk1;
        size_t zk = phi_block.zk;

        int64_t d = a;

        // optimize starting d
        if (X / (pb1 * zk1) < IACBRTX)
            d = min(d, PRIME_COUNT[X / (pb1 * zk1)]);
       
        for (; d > b + 1; --d)
        {
            assert(d < PRIMES.size());
            int64_t pd = PRIMES[d];
            // y is increasing as d is decreasing
            int64_t y = X / (pb1 * pd);
            if (y < zk1)
                continue;
            if (y >= zk)
                break;

            if (max(X / (pb1 * pb1), pb1) < pd && pd <= IACBRTX)
            {
                ++S2b;

            }
            else if (max(Z / pb1, (size_t)pb1) < pd && pd <= min( X / (pb1*pb1), IACBRTX))
            {
                
                S2b += PRIME_COUNT[X / (pb1 * pd)] - b + 1;

            }
            else
            {
                
                S2b += phi_block.sum_to(X / (pb1 * pd));
                phi_defer++;

            }

        }



        return S2b; // terminate
    }


    // Algorithm 3: computation of phi2(x,a)
    // Modification to original algorithm: use aux sieve [u,w] limited to
    // y's range in [zk1,zk) rather than acbrtx size
    int64_t P2_iter(
        const PhiBlock& phi_block,
        int64_t& vout,
        int64_t& phi_defer)
    {
        int64_t P2 = 0;
        int64_t v = 0; // thread local v

        // accumulate v
        // maintain aux sieve, u = tracking pb, y = x / u
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

        // ATOMIC ADD
        vout += v;

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
        vector<int64_t> S2(a+1);
        vector<int64_t> d2(a-1); // S2 decreasing d
        vector<char>  t(a-1);

        // Phi2
        int64_t P2 = a * (a-1) / 2; // starting sum

        // new variables
        int64_t v = a;

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

            // TODO: cleanup
            S2[b] = 0;
        }

        // create block endpoints
        // Dynamic block size of powers of 2 didn't make it faster (#5)
        // but I'm leaving in being able to specify zk endpoints
        
        vector<size_t> zks;
        for (size_t i = 1; i <= Z + BSIZE; i += BSIZE)
        {
            zks.push_back(i);
        }
        
        K = zks.size() - 1;

        // block_sum[k][b] = sum of Bk contents after sieving pb
        // = phi(zk - 1, b) - phi(zk1 - 1, b)
        vector<vector<int64_t>> phi_save(K+2, vector<int64_t>(a+1, 0));
        auto block_sum = phi_save;

        // deferred counts of phi_save from phi(y,b) calls
        auto S1_defer = phi_save; 
        auto S2_defer = phi_save;
        auto P2_defer = phi_save;
 
        // Main segmented sieve: blocks Bk = [z_{k-1}, z_k)
        // As a parallel block test: do k backwards!
        for (size_t k = K; k >= 1; --k)
        {
            // init new block
            size_t zk1 = zks[k-1];
            size_t zk = zks[k];

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
                    // TODO: ATOMIC ADD
                    S1 += S1_iter(b, phi_block, m[b], S1_defer[k][b]);
                }

                // S2 leaves, b in [astar, a-1)
                else if (astar <= b && b < a - 1)
                {
                    // don't save d2[b] between Bk
                    // TODO: ATOMIC ADD
                    S2[b] += S2_iter(b, phi_block, d2[b], t[b], S2_defer[k][b]);
                }

                // phi2, after sieved out first a primes
                else if (b == a)
                {
                    P2 += P2_iter(phi_block, v, P2_defer[k][b]); 
                }

            }
        }

        // Finalize P2
        P2 -= v*(v-1)/2;

        // sum up all deferred phi(y,b) bases
        //
        // phi_save(k,b) = full phi(zk-1,b) from Bk and all previous
        int64_t S2s = 0;

        for (size_t k = 1; k <= K; ++k)
        {
            for (size_t b = 2; b <= (size_t)a; ++b)
            {
                // accumulate phi(zk-1,b)
                phi_save[k][b] = phi_save[k-1][b] + block_sum[k][b];
                //cout << k << " " << b << " " << phi_save[k][b] << endl;
                S1    += phi_save[k-1][b] * S1_defer[k][b];
                S2[b] += phi_save[k-1][b] * S2_defer[k][b];
                P2    += phi_save[k-1][b] * P2_defer[k][b];
                
            }
        }

        for (int64_t b = 0; b <= a; ++b)
            S2s += S2[b];


        cout << "v = " << v << endl;

        // accumulate final results
        cout << "S1 = " << S1 << endl;
        cout << "S2 = " << S2s << endl;
        cout << "P2 = " << P2 << endl;

        return S0 + S1 + S2s + a - 1 - P2;
    }
};

// checks ft.sum_to(i) == sum v[0:i]
void check_ft_equal(const fenwick_tree& ft, const vector<bool>& v)
{
    int s = 0;
    for (size_t i = 0; i < v.size(); ++i)
    {
        s += v[i];
        assert(ft.sum_to(i) == s);
    }
}

void test_fenwick_tree()
{
    // example: fenwick tree over array [1,1,1,1,1]
    vector<bool> v1(5, 1);
    fenwick_tree ft(v1);
    check_ft_equal(ft, v1);

    v1[1] = 0;
    ft.try_decrease(1);
    check_ft_equal(ft, v1);
   
    v1[1] = 0;
    ft.try_decrease(1); // should not change
    check_ft_equal(ft, v1); 

    cout << "Fenwick tree tests passed" << endl;
}

// Test PhiBlock values without base match a reference
void check_phiyb(const PhiBlock& pb, const vector<int>& ref)
{
    for (size_t i = pb.zk1; i < pb.zk; ++i)
    {
        assert(pb.sum_to(i) == ref[i-pb.zk1]);
    }
}


void test_phi_block()
{
    vector<bool> ind(25, 1);
    PhiBlock pb(ind, 1, 51);
    // by design, phi block already has b = 1, p_b = 2 sieved out
    // sieved out evens, so remaining are 1,3,5,7,... = 1 mod 2
    const vector<int> phi11 =
    {
        1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
        11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,
        21,21,22,22,23,23,24,24,25,25
    };

    // sieved out by 3s, so remaining are 1, 5, 7, 11, ... = 1,5 mod 6
    vector<int> phi12 = 
    {
        1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,7,7,
        7,7,8,8,9,9,9,9,10,10,11,11,11,11,12,12,13,13,13,13,
        14,14,15,15,15,15,16,16,17,17
    };

    check_phiyb(pb, phi11);

    // sieve out b = 2, p_b = 3
    pb.sieve_out(3);
    check_phiyb(pb, phi12);

    // TODO: automated test
}




int main(int argc, char* argv[])
{
    // special test mode by supplying no arguments
    if (argc == 1)
    {
        test_fenwick_tree();
        test_phi_block();
        cout << "All tests passed!" << endl;
        return 0;
    }

    if (!(argc == 2 || argc == 4))
    {
        cerr << "Usage: ./primecount X [ALPHA BLOCKSIZE]\n";
        return 1;
    }

    // setup primecount tuning parameters to pass in

    // read float like 1e12 from command line (may not be exact for > 2^53)
    int64_t X = atof(argv[1]);
    int64_t alpha = max(1., pow(log10(X), 3) / 150); // empirical O(log^3 x)
    int64_t bsize = 1 << 20; // empirical block size

    if (argc == 4) // override defaults
    {
        alpha = atoi(argv[2]);
        bsize = atoi(argv[3]);
    }

    cout << "Computing for X = " << X << endl;
    cout << "Block size = " << bsize << endl;
    cout << "Alpha = " << alpha << endl;

    Primecount p(X, alpha, bsize);

    cout << p.primecount() << endl;
}
