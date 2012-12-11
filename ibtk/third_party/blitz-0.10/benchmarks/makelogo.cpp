#include <blitz/array.h>
#include <fstream>

using namespace blitz;

void makeLogo();

int main()
{
    makeLogo();
    return 0;
}

void setInitialConditions(Array<float,2>& c, Array<float,2>& P1, 
    Array<float,2>& P2, Array<float,2>& P3, int N, int M);

void snapshot(const Array<float,2>& P, const Array<float,2>& c);

void makeLogo()
{
    const int N = 300, M = 900;
    int niters = 3000;

    Array<float,2> P1, P2, P3, c;
    allocateArrays(shape(N,M), P1, P2, P3, c);
    Range I(1,N-2), J(1,M-2);

    setInitialConditions(c, P1, P2, P3, N, M);

    for (int iter=0; iter < niters; ++iter)
    {
        P3(I,J) = (2-4*c(I,J)) * P2(I,J)
          + c(I,J)*(P2(I-1,J) + P2(I+1,J) + P2(I,J-1) + P2(I,J+1))
          - P1(I,J);

        cycleArrays(P1,P2,P3);

        snapshot(P2, c);
    }

}

void setInitialConditions(Array<float,2>& c, Array<float,2>& P1,
    Array<float,2>& P2, Array<float,2>& P3, int N, int M)
{
    // Set the velocity field
    c = 0.3;

    ifstream ifs("blitz3.pgm");
    char tmpBuf[128];
    int pixel;
    ifs.getline(tmpBuf, 128);
    ifs.getline(tmpBuf, 128);
    ifs.getline(tmpBuf, 128);

    for (int pi=0; pi < 199; ++pi)
    {
        for (int pj=0; pj < 798; ++pj)
        {
            ifs >> pixel;
            if (pixel)
                c(pi+50,pj+56) = 0.02;
        }
    }

    // Initial pressure distribution: gaussian pulse
    using namespace blitz::tensor;
    int cr = N/6-1;
//    int cc = 7.0*M/8.0-1;
    float s2 = 64.0 * 9.0 / pow2(N/2.0);
    P1 = 0.0;
//    P2 = exp(-(pow2(i-cr)+pow2(j-cc)) * s2);
    P2 = exp(-(pow2(i-cr)) * s2);
    
    P3 = 0.0;
}


void snapshot(const Array<float,2>& P, const Array<float,2>& c)
{
    static int count = 0, snapshotNum = 0;
    if (++count < 50)
        return;

    count = 0;
    ++snapshotNum;
    char filename[128];
    sprintf(filename, "logo%03d.m", snapshotNum);

    ofstream ofs(filename);
    int N = P.length(firstDim);
    int M = P.length(secondDim);

    float Pmin = -0.6;
    float PScale = 1.0/1.2;
    float VScale = 1.0;

    ofs << "P" << snapshotNum << " = [ ";
    for (int i=0; i < N; ++i)
    {
        for (int j=0; j < M; ++j)
        {
            float value1 = (P(i,j)-Pmin)*PScale;
            float value2 = c(i,j)*VScale;
            int r1 = value1 * 4096;
            int r2 = value2 * 4096;
            ofs << r1 << " " << r2 << "    ";
        }
        if (i < N-1)
            ofs << ";" << endl;
    }
    ofs << "];" << endl;
}

