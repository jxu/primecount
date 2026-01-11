#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct {
    size_t len; // 0-based len
    long long* t; // 1-based tree, indexes [1:len]
} fenwick_tree;


long long fenwick_sum_to(fenwick_tree ft, size_t r) 
{
    long long s = 0;
    for (++r; r > 0; r -= r & -r)
        s += ft.t[r];
    return s;
}

// len is the length of the original 0-based array
void fenwick_add_to(fenwick_tree ft, size_t i, long long delta) 
{
    for (++i; i <= ft.len; i += i & -i)
        ft.t[i] += delta;
}

fenwick_tree fenwick_init(long long a[], size_t len) 
{
    fenwick_tree ft;
    ft.len = len;
    ft.t = calloc(len + 1, sizeof(long long));
    for (size_t i = 0; i < len; ++i)
        fenwick_add_to(ft, i, a[i]);
    return ft;
}


int main() 
{
    long long a[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    fenwick_tree ft = fenwick_init(a, 10);
    for (int i = 0; i < 10; ++i)
        printf("%lld\n", fenwick_sum_to(ft, i));

}