#include <iostream>
#include "popcnt_array.hpp"

int main()
{
    
    popcnt_array a(1 << 20);

    for (size_t i = 0; i < 100; ++i)
    {
        assert(a.sum_to(i) == i);
    }
    
    a.unset(20);
    a.refresh();

    for (size_t i = 0; i < 100; ++i)
        {
            std::cout << i << " " << a.sum_to(i) << std::endl;
            assert(a.sum_to(i) == (i <= 20) ? i : i-1);}
}
