#include "primecount.hpp"


int main(int argc, char* argv[])
{
    if (!(argc == 2 || argc == 5))
    {
        cerr << "Usage: ./primecount X [ALPHA BLOCKMIN BLOCKMAX]\n";
        return 1;
    }

    // setup primecount tuning parameters to pass in

    // read float like 1e12 from command line (may not be exact for > 2^53)
    double Xf = atof(argv[1]);
    if (Xf > (1ll << 53))
        cout << "WARNING: atof may not be exact, " <<
             "and you may need to change parameters for memory\n";
    if (Xf > 1e19)
        throw out_of_range("X too big!");

    // convert double to int
    u64 X = Xf;
    u64 alpha = max(1., pow(log10(X), 3) / 150); // empirical O(log^3 x)
    u64 blockmin = 16;
    u64 blockmax = 24;

    if (argc == 5) // override defaults
    {
        alpha = atoi(argv[2]);
        blockmin = atoi(argv[3]);
        blockmax = atoi(argv[4]);
    }

    cout << "Computing for X = " << X << endl;
    cout << "Alpha = " << alpha << endl;
    cout << "BLOCKMIN = " << blockmin << endl;
    cout << "BLOCKMAX = " << blockmax << endl;

    // main class init
    Primecount primecount(X, alpha, blockmin, blockmax);

    cout << primecount.primecount() << endl;
}
