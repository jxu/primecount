#include <cassert>
#include "fenwick_tree.hpp"

int main()
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

