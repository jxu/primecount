#include <cassert>
#include <iostream>
#include "bit_sieve.hpp"

int main()
{
    std::vector<bool> ind(200, 1);
    bit_sieve B(ind);

    assert(B.sum_range(0, 0) == 1);
    assert(B.sum_range(100, 100) == 1);
    assert(B.sum_range(3, 7) == 5);
    assert(B.sum_range(10, 90) == 81);
    assert(B.sum_range(0, 199) == 200);
    assert(B.sum == 200);

    B.try_decrease(100);
    assert(B.sum_range(10, 190) == 180);
    assert(B.sum == 199);

    B.try_decrease(100); // shouldn't change results
    assert(B.sum_range(10, 190) == 180);
    assert(B.sum == 199);

    std::cout << "Bitarray tests passed!" << std::endl;
}
