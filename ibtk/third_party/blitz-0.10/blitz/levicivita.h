#ifndef BZ_LEVICIVITA_H
#define BZ_LEVICIVITA_H

#include <blitz/prettyprint.h>
#include <blitz/shapecheck.h>
#include <blitz/numinquire.h>
#include <blitz/array/domain.h>
#include <blitz/etbase.h>
#include <blitz/tinyvec2io.cc>

BZ_NAMESPACE(blitz)

/* Defines an object which represents the 3-dimensional
   levi-civita symbol used for cross products. That way the cross
   product can be represented as an expression. */
class LeviCivita : public ETBase<LeviCivita> {
public:
    typedef char T_numtype;

  typedef opType<T_numtype>::T_optype T_optype;
  typedef asET<T_numtype>::T_wrapped T_typeprop;
  typedef unwrapET<T_typeprop>::T_unwrapped T_result;

    typedef void T_ctorArg1;
    typedef char       T_ctorArg2;    // dummy
  typedef TinyVector<int, 3> T_index;
  typedef LeviCivita T_range_result;
    static const int
        numArrayOperands = 0,
        numTVOperands = 0,
        numTMOperands = 0,
        numIndexPlaceholders = 3,
      minWidth = simdTypes<T_numtype>::vecWidth,
      maxWidth = simdTypes<T_numtype>::vecWidth,
        rank_ = 3;

  // copy constructor and default constructors are implicit

    int ascending(const int) const { return true; }
    int ordering(const int)  const { return INT_MIN; }
    int lbound(const int)    const { return 0; }
    int ubound(const int)    const { return 2; }

  RectDomain<rank_> domain() const
  {
    return RectDomain<rank_>(TinyVector<int,3>(0),
			     TinyVector<int,3>(2));
  }

  T_numtype operator*()   const { BZPRECONDITION(0); return 0; }
  T_numtype first_value() const { BZPRECONDITION(0); return 0; }

  bool assertInRange(const T_index& BZ_DEBUG_PARAM(index)) const {
    BZPRECHECK((index[0]<3) && (index[0]>=0) &&
	       (index[1]<3) && (index[1]>=0) &&
	       (index[2]<3) && (index[2]>=0),
	       "LeviCivita index out of range: " << index
	       << endl << "Lower bounds: 0" << endl
	       <<         "Length:       3" << endl);
    return true;
  }

#ifdef BZ_ARRAY_EXPR_PASS_INDEX_BY_VALUE
    template<int N_rank>
    T_numtype operator()(const TinyVector<int,N_rank> i) const
  { assertInRange(i); return (i[1]-i[0])*(i[2]-i[0])*(i[2]-i[1])/2; }
#else
    template<int N_rank>
    T_numtype operator()(const TinyVector<int,N_rank>& i) const
  { assertInRange(i); return (i[1]-i[0])*(i[2]-i[0])*(i[2]-i[1])/2; }
#endif

  // can't set domain
  template<int N_rank>
  const LeviCivita operator()(const RectDomain<N_rank>& d) const
  {
    BZPRECONDITION(0); return *this;
  }

    void push(int) { }
    void pop(int) { }
    void advance() { }
    void advance(int) { }
    void loadStride(int) { }

    bool isUnitStride(int) const
    { return true; }

    void advanceUnitStride()
    { }

    bool canCollapse(int,int) const
    { return true; }

    T_numtype operator[](int) const
  { BZPRECONDITION(0); return 0; }

    T_numtype fastRead(diffType) const
  { BZPRECONDITION(0); return 0; }

    T_numtype fastRead_tv(diffType) const
  { BZPRECONDITION(0); return 0; }

  // this is needed for the stencil expression fastRead to work
  void _bz_offsetData(sizeType i) const{BZPRECONDITION(0);};

    // and these are needed for stencil expression shift to work
  void _bz_offsetData(sizeType offset, int dim) const {BZPRECONDITION(0);};

  void _bz_offsetData(sizeType offset1, int dim1, sizeType offset2, int dim2) const {BZPRECONDITION(0);};

    diffType suggestStride(int) const
    { return 1; }

    bool isStride(int,diffType) const
    { return true; }

  void moveTo(int) const { BZPRECONDITION(0); }

  T_numtype shift(int offset, int dim) const {
    BZPRECONDITION(0); return T_numtype(); }

    T_numtype shift(int offset1, int dim1,int offset2, int dim2) const
  { BZPRECONDITION(0); return T_numtype();}

    template<int N_rank>
    void moveTo(const TinyVector<int,N_rank>&) const { BZPRECONDITION(0); }

    void prettyPrint(BZ_STD_SCOPE(string) &str,
        prettyPrintFormat& format) const
    {
	  str += "epsilon_ijk";
    }

    template<typename T_shape>
    bool shapeCheck(const T_shape& shape) const
  { return areShapesConformable(shape, TinyVector<int,3>(3,3,3)); }


  // we can't reduce the rank of this object, so we can't slice it
  template<typename T1, typename T2 = nilArraySection,
	   class T3 = nilArraySection, typename T4 = nilArraySection,
	   class T5 = nilArraySection, typename T6 = nilArraySection,
	   class T7 = nilArraySection, typename T8 = nilArraySection,
	   class T9 = nilArraySection, typename T10 = nilArraySection,
	   class T11 = nilArraySection>
  class SliceInfo {
  public:
    typedef void T_slice;
};

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
        typename T7, typename T8, typename T9, typename T10, typename T11>
    typename SliceInfo<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>::T_slice
    operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const
    {
      return *this;
    }
};

BZ_NAMESPACE_END

#endif
