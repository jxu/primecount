#include "primecount.hpp"

void test_fenwick_tree()
{
    // example: fenwick tree over array [1,1,1,1,1]
    fenwick_tree ft(5);

    int a1[5] = {1, 2, 3, 4, 5};

    for (int i = 0; i < 5; ++i)
        assert(ft.sum_to(i) == a1[i]);

    ft.try_decrease(1);

    int a2[5] = {1, 1, 2, 3, 4};

    for (int i = 0; i < 5; ++i)
        assert(ft.sum_to(i) == a2[i]);

    ft.try_decrease(1); // should not change

    for (int i = 0; i < 5; ++i)
        assert(ft.sum_to(i) == a2[i]);

    ft.reset();

    for (int i = 0; i < 5; ++i)
        assert(ft.sum_to(i) == a1[i]);


}

// Test PhiBlock values phi(y,b) match a reference
void test_phiyb(const PhiBlock& pb, const size_t b, const vector<int>& ref)
{
    for (size_t i = pb.zk1; i < pb.zk; ++i)
    {
        assert(pb.sum_to(i, b) == ref[i-pb.zk1]);
    }
}


void test_phi_block()
{
    const size_t bsize = 50;
    PhiBlock pb(2, bsize);
    pb.new_block(1); // block k=1: [1, 51)
    // by design, phi block already has b = 1, p_b = 2 sieved out
    // sieved out evens, so remaining are 1,3,5,7,... = 1 mod 2
    const vector<int> phi11 =
    {
        1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
        11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,
        21,21,22,22,23,23,24,24,25,25
    };

    // sieved out by 3s, so remaining are 1, 5, 7, 11, ... = 1,5 mod 6
    vector<int> phi12 = {1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,7,7,
                         7,7,8,8,9,9,9,9,10,10,11,11,11,11,12,12,13,13,13,13,
                         14,14,15,15,15,15,16,16,17,17
                        };

    int b = 1;
    test_phiyb(pb, b, phi11);
    pb.update_save(b); // prepare for next block

    // sieve out b = 2, p_b = 3
    b = 2;
    pb.sieve_out(3);
    test_phiyb(pb, b, phi12);
    pb.update_save(b); // prepare for next block

    // remaining are 51, 53, 55, ...
    vector<int> phi21 =
    {
        26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,35,
        36,36,37,37,38,38,39,39,40,40,41,41,42,42,43,43,44,44,45,45,
        46,46,47,47,48,48,49,49,50,50
    };

    // remaining 1,5 mod 6
    vector<int> phi22 =
    {
        17,17,18,18,19,19,19,19,20,20,21,21,21,21,22,22,23,23,23,23,
        24,24,25,25,25,25,26,26,27,27,27,27,28,28,29,29,29,29,30,30,
        31,31,31,31,32,32,33,33,33,33
    };

    // new block k = 2, [51,101)
    pb.new_block(2);
    b = 1;
    test_phiyb(pb, b, phi21);
    pb.update_save(b);

    b = 2;
    pb.sieve_out(3);
    test_phiyb(pb, b, phi22);
    pb.update_save(b);
}



int main()
{
    test_fenwick_tree();
    test_phi_block();
}
