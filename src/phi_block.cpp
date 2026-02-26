#include "phi_block.hpp"

// In the whole computation, Bk processed sequentially for k = 1 to K
// (potentially parallelizable by tracking phi(zk1-1,b) base values used)
// Within each Bk:
// For b from 1 to a:
//   Sieve out pb
//   Update S1b and S2b

// The physical block size is half the logical size by only storing odd values
// For example, [51, 101) would map to ind [0, 25) via y -> (y-zk1)/2

void PhiBlock::sieve_out(uint64_t pb)
{
    assert(pb > 2); // 2 already sieved by default

    // sieve out pb
    if (zk1 <= pb && pb < zk)
        phi_sum.try_decrease(tree_index(pb));

    // now only need to start at pb^2
    // (doesn't really help)
    uint64_t j0 = std::max(pb*pb, pb * ceil_div(zk1, pb));
    if (j0 % 2 == 0)
        j0 += pb; // ensure odd

    for (uint64_t j = j0; j < zk; j += 2*pb)
    {
        phi_sum.try_decrease(tree_index(j));
    }
}