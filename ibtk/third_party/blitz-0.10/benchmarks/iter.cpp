#include <iostream>
#include <blitz/array.h>
#include <blitz/timer.h>

int main() {
    using namespace blitz;

    typedef Array<double,3> Image;
    Image A(512,512,512);
    Timer timer;

    A = 0.0;
    timer.start();
    for (unsigned j=0;j<10;++j) {
        for (Image::iterator i=A.begin(),end=A.end();i!=end;++i) {
            const TinyVector<int,3> pos = i.position();
            *i += pos(0)+pos(1)+pos(2);
        }
    }
    timer.stop();
    double flops = 10.0*512*512*512*2;
    double seconds = timer.elapsedSeconds();

    double timePerOp = seconds / flops;

    cout << "ops = " << flops << endl
         << "seconds = " << seconds << endl;

    double Mflops = flops / seconds / 1.0e+6;
    cout << "Mflops = " << Mflops << endl;

    return 0;
}
