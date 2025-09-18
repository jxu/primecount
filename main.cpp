#include "primecount.hpp"

int main(int argc, char* argv[])
{
    if (!(argc == 2 || argc == 3))
    {
        cerr << "Usage: ./primecount X [ALPHA]\n";
        return 1;
    }

    // setup global constants

    // read float like 1e12 from command line (may not be exact for > 2^53)
    long X = atof(argv[1]);
    long alpha = max(1., pow(log10(X), 3) / 150); // empirical O(log^3 x)

    if (argc == 3) // override defaults
    {
        alpha = atoi(argv[2]);
    }

    // TODO: set from command line
    size_t BSIZE = 1 << 20;

    cout << "Computing for X = " << X << endl;
    cout << "Block size = " << BSIZE << endl;
    cout << "Alpha = " << alpha << endl;

    Primecount p(X, alpha, BSIZE);

    cout << p.primecount() << endl;
}
