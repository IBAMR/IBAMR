#include <blitz/timer.h>
#include <iostream>
#include <fstream>

BZ_USING_NAMESPACE(blitz)
BZ_USING_NAMESPACE(std)

extern "C" {
    void echo_f77tuned(int& N, int& niters, float& check, int& blockSize);
}

int main()
{
    int N = 1024;
    int niters = 48;
    float check;
    double Mflops = niters * 9;
    Timer timer;

    ofstream ofs("echotune.log");

    cout << "This program decides on the best block size for a typical 2D "
         << endl << "stencil operation.  Pick the block size which has the "
         << endl << "maximum Mflops/s." << endl << endl;

    cout << "Block size\tMflops/s" << endl;

    int blockSize;

    for (blockSize=1; blockSize < 32; ++blockSize)
    {
        timer.start();
        echo_f77tuned(N, niters, check, blockSize);
        timer.stop();
        cout << blockSize << "\t" << (Mflops/timer.elapsedSeconds()) << endl;
        ofs << blockSize << "\t" << (Mflops/timer.elapsedSeconds()) << endl;
    }
    for (; blockSize < 1024; blockSize += 32)
    {
        timer.start();
        echo_f77tuned(N, niters, check, blockSize);
        timer.stop();
        cout << blockSize << "\t" << (Mflops/timer.elapsedSeconds()) << endl;
        ofs << blockSize << "\t" << (Mflops/timer.elapsedSeconds()) << endl;
    }
    
    return 0;
}

