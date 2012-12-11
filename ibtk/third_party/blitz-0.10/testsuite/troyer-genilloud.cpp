// #define BZ_DEBUG_LOG_ALLOCATIONS
#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    float data[4] = { 1.0, 2.0, 3.0, 4.0 };

    {
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
        cout << endl << "Normal array ctor: " << endl;
#endif
        Array<float,1> A(4);
    }

    {
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
        cout << endl << "Old ctor version:" << endl;
#endif
        Array<float,1> A(data,shape(4));
    }

    {
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
        cout << endl << "duplicateData:" << endl;
#endif
        Array<float,1> A(data,shape(4),duplicateData);
    }

    {
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
        cout << endl << "neverDeleteData:" << endl;
#endif
        Array<float,1> A(data,shape(4),neverDeleteData);
    }

    {
#ifdef BZ_DEBUG_LOG_ALLOCATIONS
        cout << endl << "deleteDataWhenDone:" << endl;
#endif
        Array<float,1> A(new float[4], shape(4), deleteDataWhenDone);
    }

}

