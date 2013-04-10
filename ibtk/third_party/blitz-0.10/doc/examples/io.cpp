#include <blitz/array.h>
#ifdef BZ_HAVE_STD
	#include <fstream>
#else
	#include <fstream.h>
#endif

BZ_USING_NAMESPACE(blitz)

const char* filename = "io.data";

void write_arrays()
{
    ofstream ofs(filename);
    if (ofs.bad())
    {
        cerr << "Unable to write to file: " << filename << endl;
        exit(1);
    }

    Array<float,3> A(3,4,5);
    A = 111 + tensor::i + 10 * tensor::j + 100 * tensor::k;
    ofs << A << endl;

    Array<float,2> B(3,4);
    B = 11 + tensor::i + 10 * tensor::j;
    ofs << B << endl;

    Array<float,1> C(4);
    C = 1 + tensor::i;
    ofs << C << endl;
}

int main()
{
    write_arrays();

    ifstream ifs(filename);
    if (ifs.bad())
    {
        cerr << "Unable to open file: " << filename << endl;
        exit(1);
    }

    Array<float,3> A;
    Array<float,2> B;
    Array<float,1> C;

    ifs >> A >> B >> C;

    cout << "Arrays restored from file: " << A << B << C << endl;

    return 0;
}

