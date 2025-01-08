"""Computes prime_count(x) up to 10^12, along with one-time precomputation.

This module implements the algorithms of
"Computing π(x): the combinatorial method" by Oliveira e Silva.

The computation uses the Meissel method, improved by Lehmer (1959),
Lagarias-Miller-Odlyzko (1985), and Deléglise-Rivat (1996).

Basics of the Meissel-Lehmer method:
Let phi(x,a) count the positive integers up to x (inclusive) that are not
divisible by any of the first a primes, i.e. no small prime factors.
Also let phi_k(x,a) count those integers with exactly k prime factors
(including repetitions).

The Lagarias et al. methods use a = pi(alpha cbrt x), with alpha <= x^1/6
phi(x,a) = phi_0(x,a) + phi_1(x,a) + phi_2(x,a)
with phi_0(x,a) = 1, phi_1(x,a) = pi(x) - a, so that
pi(x) = phi(x,a) + a - 1 - phi_2(x,a)

The computation of phi2(x,a) is easy to do combinatorially. Since I can fit
[1,z] into memory, pi(x) can be precomputed and stored entirely for the range.

phi(x,a) is based on the recurrence which forms a binary tree:
phi(x,a) = phi(x,a-1) - phi(x/p_a,a-1)
The rest is deciding when to split the node (recurse) and when not to.
Lagarias et al. show that Lehmer's original algorithm examines too many leaves.
For efficient computation, this algorithm has lots of fancy counting of the
contributions of the leaves. Normally [1,z] needs to be sieved in segments,
but I take a shortcut as for my purposes z is small enough to fit the
whole sieve and cumulative sums (Fenwick tree for fast updates) in memory.

Of course doing it in python is a waste compared to implementing in a real
language, because I'm writing C-like code at this point.
C++ constexpr might even be able to precompute a lot at compile time.
"""

from math import comb, isqrt, prod
from itertools import accumulate


class FenwickTree:
    # Implementation is 1-based, but interface works for 0-based
    def __init__(self, l):
        self.t = [0] * (len(l)+1)
        for i in range(len(l)):
            self.add_to(i, l[i])


    def sum_to(self, r):
        """sum a[0:r]"""
        r += 1
        s = 0
        while r > 0:
            s += self.t[r]
            r -= (r & -r)
        return s


    def add_to(self, i, delta):
        """sync a[i] += delta"""
        i += 1
        while (i < len(self.t)):
            self.t[i] += delta
            i += i & -i


X_MAX = 10**6  # Module-wide maximum
ALPHA = 1       # tuning parameter, < x^1/6
C = 1           # precompute phi(x,c)

Z = (X_MAX)**(2/3) / ALPHA  # sieve max, doesn't need to be exact


def mu_pmin_sieve(N):
    """special sieve of mu(n) pmin(n) for n <= N.

    Also computes primes and PRIME_COUNT along the way.
    """
    mu_pmin = [1] * (N+1)
    primes = [0]  # 1-indexed primes
    prime_count_small = [0] * (N+1)
    pc = 0

    for j in range(2, N+1):
        if mu_pmin[j] == 1:
            for i in range(j, N+1, j):
                mu_pmin[i] = -j if mu_pmin[i] == 1 else -mu_pmin[i]

    for j in range(2, N+1):
        if mu_pmin[j] == -j:  # is prime
            primes.append(j)
            pc += 1

            for i in range(j*j, N+1, j * j):
                mu_pmin[i] = 0

        prime_count_small[j] = pc

    return mu_pmin, primes, prime_count_small


ACBRTX = int(ALPHA * X_MAX**(1/3)) + 1  # not exact
assert C <= ACBRTX

# actually more than necessary for one-off calculation pi(X_MAX)
# but benefits smaller pi(x) as in problem 501
MU_PMIN, PRIMES, PRIME_COUNT = mu_pmin_sieve(int(Z)+1)
print("sieve up to", int(Z)+1)

PRIMES_C = PRIMES[1:C+1]
Q = prod(PRIMES_C)


def pre_phi_c(N):
    """Phi_c sieves out first c primes"""
    sieve_c = [1] * (N + 1)
    sieve_c[0] = 0
    for p in PRIMES_C:
        for j in range(p, N+1, p):
            sieve_c[j] = 0

    return list(accumulate(sieve_c))


PHI_C = pre_phi_c(Q)
#print(PHI_C)


def sgn(x):
    return (x > 0) - (x < 0)


def phi_c(y):
    return (y // Q) * PHI_C[Q] + PHI_C[y % Q]


def exact_floor(z, bound_check):
    """
    Exact floor subroutine.

    :param z: initial guess (float within 1 of the true result)
    :param bound_check: exact predicate of z satisfying bound
    :return: largest integer satisfying bound

    https://cs.stackexchange.com/a/4841
    """
    z = int(z) + 1
    while not bound_check(z):
        z -= 1
    return z


def prime_count(x):
    if x < 0:
        raise ValueError

    if x <= Z:
        return PRIME_COUNT[x]

    assert ALPHA <= x**(1/6)
    acbrtx = ALPHA * x**(1/3)  # float, may be off
    z = int(x / acbrtx) + 1  # not exact

    # find exact floor of acbrtx
    iacbrtx = exact_floor(acbrtx, lambda z: z**3 <= ALPHA**3 * x)
    
    print(iacbrtx, isqrt(x))

    a = PRIME_COUNT[iacbrtx]  # pi(alpha x^1/3)
    a2 = PRIME_COUNT[isqrt(x)]  # pi(x^1/2)

    print(f"a {a} a2 {a2}")

    # phi2(x,a)
    P2 = comb(a, 2) - comb(a2, 2)
    for b in range(a+1, a2+1):
        # if PRIME_COUNT up to z not available, use phi(x/p_b,a) + a - 1
        # from phi sieved below up to a. Each query takes log(z) time
        P2 += PRIME_COUNT[x // PRIMES[b]]


    # phi(x,a) Ordinary leaves, only alpha cbrt x of them
    S0 = 0
    for n in range(1, iacbrtx+1):
        # pmin(1) = +inf
        if n == 1 or abs(MU_PMIN[n]) > PRIMES[C]:
            S0 += sgn(MU_PMIN[n]) * phi_c(x // n)


    # phi(x,a) Special leaves
    S = 0
    # sieve out first C primes, only storing odd values
    assert C >= 1
    sieve_ind = [1] * (z // 2)

    for p in PRIMES_C[1:]:  # skip 2
        for j in range(p, z, 2*p):  # odd multiples of p
            sieve_ind[j//2] = 0

    # phi(y,b) = tree sum up to index (y-1)//2
    dyn_sieve = FenwickTree(sieve_ind)

    S1_total = 0

    for b in range(C, a-1):
        #print(b)
        pb1 = PRIMES[b+1]

        # Case 1 leaves, Algorithm 1
        if pb1**2 <= iacbrtx:
            #print("case1")
            S1b = 0
            m1b = iacbrtx
            while (m1b * pb1)**3 > ALPHA**3 * x:
                if abs(MU_PMIN[m1b]) > pb1:
                    y = x // (m1b * pb1)
                    phi_b = dyn_sieve.sum_to((y-1)//2)
                    S1b -= sgn(MU_PMIN[m1b]) * phi_b
                m1b -= 1

            S += S1b
            S1_total += S1b

        # Case 2 leaves, Algorithm 2 hell
        else:
            xpb12 = x // (pb1**2)

            # number of trivial leaves is a + 1 - tb
            if xpb12 <= pb1:        tb = b + 2
            elif xpb12 < iacbrtx:   tb = PRIME_COUNT[xpb12] + 1
            else:                   tb = a + 1

            # step 1
            d2b = tb - 1  # largest d not considered yet
            S2b = a - d2b
            t = 0

            while d2b > b + 1:  # step 2
                y = x // (pb1 * PRIMES[d2b])

                if t == 0:  # step 3, clustered easy leaves
                    if y >= iacbrtx:
                        t = 2
                    else:
                        l = PRIME_COUNT[y] - b + 1
                        d_ = PRIME_COUNT[x // (pb1 * PRIMES[b + l])]

                        # step 4
                        if PRIMES[d_+1]**2 <= x // pb1 or d_ <= b:
                            t = 1
                            # goto step 6
                        else:
                            S2b += l * (d2b - d_)
                            d2b = d_
                            continue  # goto 2

                if t == 1:  # step 5, sparse easy leaves
                    if y >= iacbrtx:
                        t = 2
                    else:
                        l = PRIME_COUNT[y] - b + 1
                        # step 6
                        S2b += l
                        d2b -= 1
                        continue  # goto 2

                if t == 2:  # step 7-9, hard leaves
                    S2b += dyn_sieve.sum_to((y-1)//2)
                    d2b -= 1

            S += S2b


        # sieve out (odd) multiples of p_(b+1) for next step
        for i in range(pb1, z, 2*pb1):
            if sieve_ind[i//2] == 1:
                dyn_sieve.add_to(i//2, -1)
                sieve_ind[i//2] = 0

    print("S1total", S1_total)
    print(a, P2, S0, S) 
    return S0 + S + a - P2 - 1


def test_prime_count():
    # A006880
    powers_10 = (0, 4, 25, 168, 1229, 9592, 78498, 664579, 5761455, 50847534,
                 455052511, 4118054813, 37607912018)

    for i in range(len(powers_10)):
        assert prime_count(10**i) == powers_10[i]

print(prime_count(X_MAX))
