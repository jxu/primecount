#include "primecount.hpp"
#include <stdexcept>


int main(int argc, char* argv[])
{
    if (!(argc == 2 || argc == 5))
    {
        std::cerr << "Usage: ./primecount X [ALPHA BLOCKMIN BLOCKMAX]\n";
        return 1;
    }

    // setup primecount tuning parameters to pass in

    // read float like 1e12 from command line (may not be exact for > 2^53)
    double Xf = atof(argv[1]);
    if (Xf > (1ll << 53))
        std::cout << "WARNING: atof may not be exact, " <<
            "and you may need to change parameters for memory\n";
    if (Xf > 1e19)
        throw std::out_of_range("X too big!");

    // convert double to int
    uint64_t X = Xf;
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
