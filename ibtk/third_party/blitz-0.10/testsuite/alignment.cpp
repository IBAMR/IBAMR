#include "testsuite.h"
#include <blitz/array.h>
#include <vector>

BZ_USING_NAMESPACE(blitz)

template<typename T_numtype, int N, int M>
void testalign() {
  BZTEST(__alignof__(TinyVector<T_numtype,N>)==simdTypes<T_numtype>::byteWidth);
  BZTEST(__alignof__(TinyMatrix<T_numtype,N,M>)==simdTypes<T_numtype>::byteWidth);

  const int w=simdTypes<T_numtype>::vecWidth;
  BZTEST(sizeof(TinyVector<T_numtype,N>)==(N/w + (N%w?1:0))*w*sizeof(T_numtype));
  BZTEST(sizeof(TinyMatrix<T_numtype,N,M>)==(N*M/w + (N*M%w?1:0))*w*sizeof(T_numtype));
}

int main()
{
// test that TinyVector and TinyMatrix are aligned properly.
  testalign<int,1,1>();
  testalign<int,2,2>();
  testalign<int,3,2>();
  testalign<int,4,2>();
  testalign<int,5,2>();
  testalign<int,6,2>();
  testalign<int,7,2>();

  testalign<float,1,2>();
  testalign<float,2,2>();
  testalign<float,3,2>();
  testalign<float,4,2>();
  testalign<float,5,2>();
  testalign<float,6,2>();
  testalign<float,7,2>();

  testalign<double,1,2>();
  testalign<double,2,2>();
  testalign<double,3,2>();
  testalign<double,4,2>();
  testalign<double,5,2>();
  testalign<double,6,2>();
  testalign<double,7,2>();

  // test that simdTypes returns correct values

  BZTEST(bzCC(_bz_meta_bitwidth<0,0>::width==0));
  BZTEST(bzCC(_bz_meta_bitwidth<1,0>::width==1));
  BZTEST(bzCC(_bz_meta_bitwidth<2,0>::width==2));
  BZTEST(bzCC(_bz_meta_bitwidth<3,0>::width==2));
  BZTEST(bzCC(_bz_meta_bitwidth<4,0>::width==3));
  BZTEST(bzCC(_bz_meta_bitwidth<5,0>::width==3));
  BZTEST(bzCC(_bz_meta_bitwidth<6,0>::width==3));
  BZTEST(bzCC(_bz_meta_bitwidth<7,0>::width==3));
  BZTEST(bzCC(_bz_meta_bitwidth<8,0>::width==4));

  for(size_t i=0; i<100; ++i) {
    const size_t l = simdTypes<float>::paddedLength(i);
    BZTEST(l%simdTypes<float>::vecWidth==0);
    BZTEST(l>=i);
    BZTEST(l-i<simdTypes<float>::vecWidth);
  }

  TinyVector<float, 2*simdTypes<float>::vecWidth> v;
  const int w=simdTypes<float>::vecWidth;

  BZTEST(simdTypes<float>::offsetToAlignment(v.data())==0);
  for (int i=1; i<w; ++i) {
    BZTEST(simdTypes<float>::offsetToAlignment(v.data()+i)==w-i);
  }

  if(w>1) {
    // test that arrays correctly report as being aligned or not.
    Array<float, 2> A(3*w,3*w);
    BZTEST(A.isVectorAligned(0));
    Array<float,1> B(A(1,Range::all()));
    BZTEST(B.isVectorAligned(0));
    BZTEST(!A(Range::all(),1).isVectorAligned(0));
    BZTEST(B(Range(w,2*w-1)).isVectorAligned(0));
    BZTEST(!B(Range(1,2*w-1)).isVectorAligned(0));
    BZTEST(B(Range(1,2*w-1)).isVectorAligned(w-1));

    // test that expressions correctly report as being aligned or not
    BZTEST((A(1,Range::all())*B).isVectorAligned(0));
    BZTEST(!(A(Range::all(),1)*B).isVectorAligned(0));
    BZTEST((A(Range::all(),1)*B(Range(1,2*w-1))).isVectorAligned(w-1));

    // test that we don't overwrite the end of aligned but uneven-length arrays
    B=0;
    B(Range(0,0))=1;
    BZTEST(sum(B)==1);
    B(Range(0,w))=1;
    BZTEST(sum(B)==w+1);
  }

  // This mixed-type vector assignment caused an alignment fault if
  // simd width was set to one, because then we did not have any
  // alignment but used aligned instructions. However, it only
  // segfaulted on O3, so not by default in the testsuite either
  {    
    std::vector<TinyVector<double,3> > a;
    a.resize(10);
    TinyVector<float,3> b;
    cout << *reinterpret_cast<long int*>(&a[0])%16 << '\t' << &b << endl;
    cout << *reinterpret_cast<long int*>(&a[1])%16 << '\t' << &b << endl;
    a[0]=b;
    a[1]=b;
  }
  
  return 0;
}

