#include "primecount.hpp"


int main(int argc, char* argv[])
{
    if (!(argc == 2 || argc == 4))
    {
        cerr << "Usage: ./primecount X [ALPHA BLOCKSIZE]\n";
        return 1;
    }

    // setup primecount tuning parameters to pass in

    // read float like 1e12 from command line (may not be exact for > 2^53)
    int64_t X = atof(argv[1]);
    int64_t alpha = max(1., pow(log10(X), 3) / 150); // empirical O(log^3 x)
    int64_t bsize = 1 << 24; // empirical block size

    if (argc == 4) // override defaults
    {
        alpha = atoi(argv[2]);
        bsize = atoi(argv[3]);
    }

    cout << "Computing for X = " << X << endl;
    cout << "Block size = " << bsize << endl;
    cout << "Alpha = " << alpha << endl;

    Primecount p(X, alpha, bsize);

    cout << p.primecount() << endl;
}
