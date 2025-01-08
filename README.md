## Basics of the Meissel-Lehmer method

This program implements the algorithms of
"Computing π(x): the combinatorial method" by Oliveira e Silva (2006).

Legendre (1808) was the first to notice $\pi(x)$ does not require explicitly determining all primes up to $x$, making use of inclusion-exclusion. The computation uses the Meissel (1870) method, improved by Lehmer (1959),
Lagarias-Miller-Odlyzko (1985), and Deléglise-Rivat (1996). 

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

The algorithm as a whole takes $O(x^{2/3})$ time, and the sieve is segmented to only use $O(\alpha \sqrt[3] x)$ memory. Without segmenting, it is more straightforward to sieve the entire interval $[1,x^{2/3} / \alpha)$. Over each interval $B_k$, the contribution of $\phi(y,b)$ for $y \in B_k$ and $c \le b < a$ is computed.

