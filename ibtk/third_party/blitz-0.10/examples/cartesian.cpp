#include <blitz/array.h>
#include <blitz/array/cartesian.h>
#include <list>

BZ_USING_NAMESPACE(blitz)

int main()
{
    list<int> index1;
    index1.push_back(5);
    index1.push_back(7);
    index1.push_back(9);

    list<int> index2;
    index2.push_back(1);
    index2.push_back(2);
    index2.push_back(4);
    index2.push_back(7);

    CartesianProduct<TinyVector<int,2>,list<int>,2> indexSet(index1,index2);
    typedef CartesianProduct<TinyVector<int,2>,list<int>,2>::iterator T_iterator;

    T_iterator iter = indexSet.begin(), end = indexSet.end();
    for (; iter != end; ++iter)
        cout << *iter << endl;

    return 0;
}

