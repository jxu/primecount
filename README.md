# Primecount

This program implements the algorithms of
["Computing π(x): the combinatorial method" by Tomás Oliveira e Silva (2006)](https://sweet.ua.pt/tos/bib/5.4.pdf). 
It is also a learning project for numerical computing in C++ and a sandbox to test out optimizations.

## Tuning and Benchmarks

Alpha is a tuning parameter that trades fewer segmented sieve blocks to iterate over for more computation per block. From the paper, it must be <= x^(1/6) and should grow like O(log^3 x). I found log10(x)^3 / 150 generated reasonable alpha values (in the single-threaded original implementation), and there is a decent amount of leeway.

The following benchmarks used a block size of 2^16 to 2^24 and default alpha value. 
Compiled with g++ -O3 and run on my laptop (Ubuntu 24.04, GCC 13.3.0, i7-7700HQ with 4 cores, 24 GB memory).
1 thread was run with `OMP_NUM_THREADS=1`.

| x     | Time (1 thread) | Time (8 threads) | π(x)               |
|-------|-----------------|------------------|--------------------|
| 10^12 | 0.5s            | 0.3s             | 37607912018        |
| 10^13 | 1.4s            | 0.6s             | 346065536839       |
| 10^14 | 5.6s            | 2.0s             | 3204941750802      |
| 10^15 | 22s             | 7.7s             | 29844570422669     |
| 10^16 | 1m36s           | 31s              | 279238341033925    |
| 10^17 | 6m52s           | 2m17s            | 2623557157654233   |
| 10^18 | 29m             | 9m48s            | 24739954287740860  |
| 10^19 | 186m*           | 45m              | 234057667276344607 |

*Computed on a non-parallel older version which fit into memory.

*Still a better use of electricity than crypto*

## Correctness

The powers of 10 in the table match Table IV of the paper.
I use integer math mostly, but the most likely errors if any will be from imprecise floating point calculations like cbrt(x), where I'm not sure if the exact floor integer is needed.
Previously I had some functions to calculate these exactly by checking with integer math, but I removed those.
Maybe this causes some rare off-by-one error somewhere.

Since I capped the max X input value, overflow shouldn't be an issue.
Incidentally, 10^19 is just outside of the max of a signed 64-bit int (approx. 9.22 * 10^18), so I have a mess of unsigned and signed ints everywhere.

## Optimizations

- Odd phi block (mentioned in Algorithm 4): simple to index and saves block space to fit into cache
- ~Fenwick tree sign bit instead of bool array (Appendix):~ saves a little space, but not faster
- ~Mod 3, 5 wheel in phi block (my idea):~ saves a little space, but not faster and adds a lot of indexing complexity
- Dynamic block size (section II-D): This was actually slower for single-threaded, but in multi-threading it helps a little to balance out the block work for smaller X
- Phi block initialization from phi(n,c) table (Algorithm 4): not actually faster but I'm keeping it
- Parallel blocks (not covered in paper): see below

Need more tests!

## Parallelization

After my other optimizations failed, the remaining obvious optimization was multi-threading, with the alluring promise of easy large speedups (although it comes with its own challenges, like false sharing). 

The way the algorithms are written in the paper assume sequentially going through the blocks from block 1 to block K.
Only the starting variables depend on the last block as they pick up where the last left off.

The simple way around this I thought of was to make each block completely independent by using a fixed starting point, computed from the bounds the block endpoints imply.
The block is only used for phi(y,b), so y being in the block is the main constraint and checked when needed.
Now, the value of phi(y,b) depends on this block plus summing the phi values from previous blocks, but because I don't have those previous numbers, I defer those additions and keep track of what would've been added until the end, when they are added in.

* For Algorithm 3 (computation of phi2), the author uses an auxiliary sieve to sieve portions of iacbrtx size "on-demand", that is sieve a new portion when we've exhausted the current portion.
Instead, I chose to sieve the range of `u` that corresponded to the `y` in the block.
The auxiliary sieve in total only goes up to $\sqrt x$, so it is comparatively less work than S2.

* For Algorithm 2 (computation of S2b), the goto structure as presented is confusing, so I did the leaves in a straightforward manner from the three types (trivial, easy, hard).
I was not able to incorporate the optimization for the clustered easy leaves because of the jump in setting $d_{2b} = d'$.
In a sequential computation, the exact $d_{2b}$ is set for the next block, while in mine $d_{2b}$ as a starting point is estimated from some of the constraints.
I think if $d'$ jumps somewhere in the next block, there is no way of accounting for what it skipped over, while going through $d$ one-by-one ensures no $d$ is missed and exactly the block interval is considered.
From some quick testing it looked like the clustered easy leaves optimization only saved a little bit of time anyway, like 10%.

* For X = 10^19, the memory usage of the phi and S2 defer tables was too large, as each take about a * K * 8 bytes.
While a can be reduced by alpha and K can be reduced by a larger block size, it still took too much memory.
My workaround was to not parallelize k completely, but only parallelize a batch of k into say L = 64 blocks, then regroup and collect sums, and repeat.
Now each table takes only a * L * 8 bytes, where a = O(acbrtx).
The whole program relies on having a full prime count table up to acbrtx, so O(acbrtx) is the storage required without overhauling the entire math.


## Basics of the Meissel-Lehmer method

Legendre (1808) was the first to notice $\pi(x)$ does not require explicitly determining all primes up to $x$, making use of inclusion-exclusion. The computation uses the Meissel (1870) method, improved by Lehmer (1959),
Lagarias-Miller-Odlyzko (1985), and Deléglise-Rivat (1996). 
It is combinatorial in nature and requires only elementary number theory, with some basics about prime factorization and the Möbius function's role in inclusion-exclusion. 

Lagarias and Odlyzko (1987) described an analytic algorithm with better time complexity O(x^(1/2+ε)), but it goes way over my head and in practice has not been faster than the combinatorial method yet. [The implementation](https://arxiv.org/pdf/1203.5712) also uses interval arithmetic. Don't expect that project from me any time soon.

Let $\phi(x,a)$ count the positive integers up to $x$ (inclusive) that are not divisible by any of the first $a$ primes, i.e. no small prime factors. Also let $\phi_k(x,a)$ count those integers with exactly $k$ prime factors (including repetitions).

By the fundamental theorem of arithmetic, 

$$\phi(x,a) = \phi_0(x,a) + \phi_1(x,a) + \phi_2(x,a) + \cdots$$

Since the prime factors have to be greater than $p_a$, we have $\phi_k(x,a)=0$ when $a \ge \pi (\sqrt[k] x)$. We also have $\phi_0(x,a) = 1$ (the only number with $0$ prime factors is $1$), and for $a \le \pi(x)$ we have $\phi_1(x,a) = \pi(x) - a$ (the primes $> p_a$). 

The Lagarias et al. methods use $a = \pi(\alpha \sqrt[3]x)$, with $\alpha \le x^{1/6}$, thus $a \le \pi(\sqrt x)$. Therefore 

$$\phi(x,a) = \phi_0(x,a) + \phi_1(x,a) + \phi_2(x,a)$$

rearranging

$$\pi(x) = \phi(x,a) + a - 1 - \phi_2(x,a)$$

The computation of $\phi_2(x,a)$ is relatively easy to do combinatorially. 

The main computation $\phi(x,a)$ is based on the recurrence which forms a binary tree:

$$
\begin{align*}
\phi(x,a) &= \sum_{n=1}^{\lfloor x \rfloor} ([p_{\min}(n)  \ge p_a] - [p_{\min}(n) = p_a]) \\
&= \phi(x,a-1) - \sum_{n=1}^{\lfloor x / p_a \rfloor} [p_{\min}(n) \ge p_a] \\
&= \phi(x,a-1) - \phi(x/p_a, a-1)
\end{align*} 
$$

Lehmer's original algorithm stops recursing at $\phi(x,c)$ for a pre-computed table with fixed $c$. Lagarias et al. show that  examines too many leaves, so the rest of the improvements are deciding when to split the node (recurse) and when not to (leaf), along with lots of fancy counting of the contributions of the leaves using a sieve. 

The algorithm as a whole takes $O(x^{2/3})$ time, and the sieve is segmented to only use $O(\alpha \sqrt[3] x)$ memory. Without segmenting, it is more straightforward to sieve the entire interval $[1,z)$ where $z = x^{2/3} / \alpha$. Over each interval $B_k = [z_{k-1},z_k)$, the contribution of $\phi(y,b)$ for $y \in B_k$ and $c \le b < a$ is computed. A Fenwick tree is used to be able to compute this as a prefix sum (plus saved $\phi(z_{k-1}-1,b)$ from previous blocks) efficiently while also allowing fast updates to an indicator array from sieving out $p_b$.

