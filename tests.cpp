#include "primecount.hpp"

using namespace std;

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
    fenwick_tree ft(5);
    check_ft_equal(ft, v1);

    v1[1] = 0;
    ft.try_decrease(1);
    check_ft_equal(ft, v1);
   
    v1[1] = 0;
    ft.try_decrease(1); // should not change
    check_ft_equal(ft, v1); 

    v1.assign(5, 1); // reset 
    ft.reset();
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
    PhiBlock pb(1, 51);
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



int main()
{
    test_fenwick_tree();
    test_phi_block();
}
