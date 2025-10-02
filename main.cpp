#include "primecount.hpp"


int main(int argc, char* argv[])
{
    if (!(argc == 2 || argc == 5))
    {
        std::cerr << "Usage: ./primecount X [ALPHA BLOCKMIN BLOCKMAX]\n";
        return 1;
    }

    // setup primecount tuning parameters to pass in

    // read float like 1e12 from command line (may not be exact for > 2^53)
    int64_t X = atof(argv[1]);
    int64_t alpha = std::max(1., pow(log10(X), 3) / 150); // empirical O(log^3 x)
    int64_t blockmin = 16;
    int64_t blockmax = 24;

    if (argc == 5) // override defaults
    {
        alpha = atoi(argv[2]);
        blockmin = atoi(argv[3]);
        blockmax = atoi(argv[4]);
    }

    std::cout << "Computing for X = " << X << std::endl;
    std::cout << "Alpha = " << alpha << std::endl;
    std::cout << "BLOCKMIN = " << blockmin << std::endl;
    std::cout << "BLOCKMAX = " << blockmax << std::endl;

    // main class init
    Primecount primecount(X, alpha, blockmin, blockmax);

    std::cout << primecount.primecount() << std::endl;
}
