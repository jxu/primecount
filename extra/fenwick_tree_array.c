#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

typedef long long* fenwick_tree; // plain array data structure

long long fenwick_sum_to(fenwick_tree ft, size_t r) 
{
    long long s = 0;
    for (++r; r > 0; r -= r & -r)
        s += ft[r];
    return s;
}

// len is the length of the original 0-based array
void fenwick_add_to(fenwick_tree ft, size_t len, size_t i, long long delta) 
{
    for (++i; i <= len; i += i & -i)
        ft[i] += delta;
}

// array can init globally if size is known 
fenwick_tree fenwick_init(long long a[], size_t len) 
{
    fenwick_tree ft = calloc(len + 1, sizeof(long long));
    for (size_t i = 0; i < len; ++i) 
        fenwick_add_to(ft, len, i, a[i]);
    return ft;
}

int main() 
{
    long long a[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    fenwick_tree ft = fenwick_init(a, 10);

    for (int i = 0; i < 10; ++i)
        printf("%lld\n", fenwick_sum_to(ft, i));
}