#include "phi_block.h"
#include <stdio.h>

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

    std::cout << "Fenwick tree tests passed" << std::endl;
    */
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

int main()
{
    test_fenwick_tree();
    test_phi_block();
    printf("All tests passed!\n");
}

