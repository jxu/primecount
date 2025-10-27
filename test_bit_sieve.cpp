#include <cassert>
#include <iostream>
#include <random>
#include "bit_sieve.hpp"

int main()
{
    // test all 0s
    std::vector<bool> ind0(200, 0);
    bit_sieve B0(ind0);

    assert(B0.sum_range(0, 5) == 0);

    // test all 1s and decreasing index
    std::vector<bool> ind(200, 1);
    bit_sieve B(ind);

    assert(B.sum_range(0, 0) == 1);
    assert(B.sum_range(0, 63) == 64); // better test for rmask shifting
    assert(B.sum_range(100, 100) == 1);
    assert(B.sum_range(3, 7) == 5);
    assert(B.sum_range(10, 90) == 81);
    assert(B.sum_range(0, 199) == 200);
    assert(B.sum == 200);

    B.try_decrease(100);
    assert(B.sum_range(10, 190) == 180);
    assert(B.sum_range(10, 20) == 11);
    assert(B.sum == 199);

    B.try_decrease(100); // shouldn't change results
    assert(B.sum_range(10, 190) == 180);
    assert(B.sum_range(10, 20) == 11);
    assert(B.sum == 199);

    // randomized testing
    std::mt19937 rng(1229);
    const int n = 256;
    std::uniform_int_distribution<> unif1(0, 1);
    std::uniform_int_distribution<> unifn(0, n-1);

    for (int t = 0; t < 10; ++t) // trials
    {
        // randomly fill bool vector
        std::vector<bool> ind(n);
        for (int j = 0; j < n; ++j)
            ind[j] = unif1(rng);

        // init tree from vector
        bit_sieve B(ind);

        // pick random indices to decrease and test sums
        for (int j = 0; j < 100; ++j)
        {
            int x = unifn(rng);
            ind[x] = 0;
            B.try_decrease(x);

            for (int k = 0; k < 100; ++k)
            {
                int l = unifn(rng);
                int r = unifn(rng);

                if (r < l)
                    std::swap(l, r);

                int s = 0;
                for (int i = l; i <= r; ++i)
                    s += ind[i];

                assert(B.sum_range(l, r) == s);

            }
        }
    }
    

    std::cout << "Bitarray tests passed!" << std::endl;
}
