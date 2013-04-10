#!/bin/env python

# Generates stencil code. This replaces the macros in stencil-et.h,
# which make it impossible to debug the generated code.

import sys

os=open(sys.argv[1],'w')
print "Generating file %s"%sys.argv[1]

def BZ_ET_STENCIL_REDIRECT(name):
    stub="""


/* Explicit operators for arrays for stencil name. */
template<typename T, int N>						
inline _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> >
#name#(const Array<T,N>& d1)
{ return #name#(d1.wrap()); }

template<typename T, int N>
inline _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> >
#name#(Array<T,N>& d1)
{ return #name#(d1.wrap()); }
"""
    os.write(stub.replace("#name#",name))

def BZ_ET_STENCIL(name, result, etresult, MINB, MAXB):
    stub="""


/** Defines a stencil ET "#name#" that operates on an Array<P_numtype, N_rank>
   and specifies the return type as Array<#result#, N_rank>. The result
   type is used when running on an array and the etresult type when
   running on an expression. The extent of the stencil is MINB-MAXB.
   If you want to refer to the native type
   of the expression, set result="P_numtype" and etresult="typename
   T1::T_numtype". Sorry for that ugliness, but they define types
   differently. The stencil ET calls the stencil operator
   name_stencilop, defined in stencilops.h. **/

  template<typename P_expr, _bz_typename P_numtype>			
  class #name#_et : public _bz_StencilExpr<P_expr, P_numtype>		
  {									
  public:								
    typedef _bz_StencilExpr<P_expr, P_numtype> T_base;			
    typedef _bz_typename T_base::T_numtype T_numtype;			
    typedef _bz_typename T_base::T_expr T_expr;				

  // if P_numtype is an ET-type, we need to return an expr
  typedef typename selectET<P_numtype,
			    T_numtype, 
      ETBase<_bz_ArrayExpr<_bz_ArrayExprConstant<P_numtype> > > >::T_selected T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;
  typedef T_numtype T_optype;

  template<int N> struct tvresult {
    typedef #name#_et<
      typename T_expr::template tvresult<N>::Type,
      T_numtype> Type;
  };

    typedef  #name#_et<_bz_typename P_expr::T_range_result, T_numtype> T_range_result; 
									
    using T_base::iter_;						
    using T_base::rank_;							
  public:								
    #name#_et(const #name#_et& a) :					
    _bz_StencilExpr<P_expr, T_numtype>(a)				
      { }								
									
    #name#_et(BZ_ETPARM(T_expr) a) :					
    _bz_StencilExpr<P_expr, T_numtype>(a)				
      { }								
									
    #name#_et(_bz_typename T_expr::T_ctorArg1 a) :			
    _bz_StencilExpr<P_expr, T_numtype>(a)				
      { }								
									
    T_result operator*() const						
    { return #name#_stencilop(iter_); }						
									
    /* this is not really const, because we don't undo the moveTo, but	
       that should not be visible to outside. It would be if we used	
       some kind of mixed index and stack traversal, but that will	
       screw things up, const or not. */				
    template<int N_rank2>						
      T_result operator()(const TinyVector<int, N_rank2>& i) const	
    { iter_.moveTo(i); return #name#_stencilop(iter_); }    									
    T_range_result operator()(const RectDomain<rank_>& d) const		
    { return T_range_result(iter_(d)); }				
									
    T_result operator[](int i) const					
    { return #name#_stencilop(iter_[i]); }									
    T_result fastRead(diffType i) const				
    {/* this probably isn't very fast... */				
      iter_._bz_offsetData(i);						
      T_result r = #name#_stencilop(iter_);					
      iter_._bz_offsetData(-i);						
      return r;								
    }									

    /** This way of vectorizing won't work on stencils. */
    template<int N>
    typename tvresult<N>::Type fastRead_tv(diffType i) const {
      BZPRECHECK(0, "Can't vectorize stencils");
      return iter_.fastRead_tv<N>(i); }
      
    T_result shift(int offset, int dim) const				
    {									
      iter_._bz_offsetData(offset, dim);				
      T_result r = #name#_stencilop(iter_);					
      iter_._bz_offsetData(-offset, dim);				
      return r;								
    }									
									
    T_result shift(int offset1, int dim1, int offset2, int dim2) const	
    {									
      iter_._bz_offsetData(offset1, dim1, offset2, dim2);		
      T_result r = #name#_stencilop (iter_);					
      iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		
      return r;								
    }									
									
    void prettyPrint(BZ_STD_SCOPE(string) &str,				
		     prettyPrintFormat& format) const			
    {									
      str += "#name# (stencil)";						
      str += "(";							
      iter_.prettyPrint(str, format);					
      str += ")";							
    }									
    									
    template<typename T1, typename T2 = nilArraySection,		
      class T3 = nilArraySection, typename T4 = nilArraySection,	
      class T5 = nilArraySection, typename T6 = nilArraySection,	
      class T7 = nilArraySection, typename T8 = nilArraySection,	
      class T9 = nilArraySection, typename T10 = nilArraySection,	
      class T11 = nilArraySection>					
      class SliceInfo {							
      public:								
      typedef #name#_et<T_expr, T_numtype> T_slice;			
      };								
    									
    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
      typename T7, typename T8, typename T9, typename T10, typename T11> 
      #name#_et							
      operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const 
    {									
    /* because stencils work inherently in several dimensions it's	
       complicated to slice the domain. slices will be changed to unit	
        ranges instead. slicing stencil result thus *never* changes	
       the rank of the expression, unlike the normal case. */		
      return #name#_et						
	(iter_(_bz_makeRange(r1),					
	       _bz_makeRange(r2),					
	       _bz_makeRange(r3),					
	       _bz_makeRange(r4),					
	       _bz_makeRange(r5),					
	       _bz_makeRange(r6),					
	       _bz_makeRange(r7),					
	       _bz_makeRange(r8),					
	       _bz_makeRange(r9),					
	       _bz_makeRange(r10),					
	       _bz_makeRange(r11)));					
    }									
									
  };									
  /* generate an ET object from an expression */			
  template<typename T1>							
  inline _bz_ArrayExpr<#name#_et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, #etresult#> > 
  #name#(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1)				
  {									
    return _bz_ArrayExpr<#name#_et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, #etresult#> > 
      (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),MINB, MAXB))); 
  }									
  /* redirect calls with bare arrays to the main function */		
  template<typename T, int N>						
  inline _bz_ArrayExpr<#name#_et<typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, #result#> > 
  #name#(const Array<T,N>& d1)						
  { return #name#(d1.wrap()); }						
									
  template<typename T, int N>						
  inline _bz_ArrayExpr<#name#_et<typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, #result#> > 
  #name#(Array<T,N>& d1)							
   { return #name#(d1.wrap()); }"""
    os.write(stub.replace("#name#",name).replace("#result#",result).replace("#etresult#",etresult).replace("MINB",str(MINB)).replace("MAXB",str(MAXB)))


def BZ_ET_STENCIL2(name, result, etresult, MINB, MAXB):
    stub="""


/** Defines a stencil "#name#" ET that operates on two arrays of arbitrary type
   and specifies the return type as array<#result#, N_rank>. The result
   type is used when running on an array and the etresult type when
   running on an expression. If you want to refer to the native type
   of the expression, set result="P_numtype" and etresult="typename
   T1::T_numtype". Sorry for that ugliness, but they define types
   differently. */
template<typename P_expr1, typename P_expr2, _bz_typename P_numtype>	
class #name#_et2 : public _bz_StencilExpr2<P_expr1, P_expr2, P_numtype> 
{									
public:								
  typedef _bz_StencilExpr2<P_expr1, P_expr2, P_numtype> T_base;	
  typedef _bz_typename T_base::T_numtype T_numtype;			
  typedef _bz_typename T_base::T_expr1 T_expr1;			
  typedef _bz_typename T_base::T_expr2 T_expr2;			

  // if P_numtype is an ET-type, we need to return an expr
  typedef typename selectET<P_numtype,
			    T_numtype, 
      ETBase<_bz_ArrayExpr<_bz_ArrayExprConstant<P_numtype> > > >::T_selected T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;
  typedef T_numtype T_optype;

  /// dummy
  template<int N> struct tvresult {
    typedef #name#_et2<
      typename T_expr1::template tvresult<N>::Type,
      typename T_expr2::template tvresult<N>::Type,
      T_numtype> Type; 
  };

  typedef  #name#_et2<_bz_typename P_expr1::T_range_result, _bz_typename P_expr2::T_range_result, T_numtype> T_range_result; 
									
  using T_base::iter1_;						
  using T_base::iter2_;						
  using T_base::rank_;							
public:								
  #name#_et2(const #name#_et2& a) :				
  _bz_StencilExpr2<P_expr1, P_expr2, T_numtype>(a)			
  { }								
									
  #name#_et2(BZ_ETPARM(T_expr1) a, BZ_ETPARM(T_expr2) b) :		
  _bz_StencilExpr2<P_expr1, P_expr2, T_numtype>(a, b)		
  { }								

  T_result operator*() const						
  { return #name#_stencilop(iter1_, iter2_); }					
									
  T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const 
  { iter1_.moveTo(i); iter2_.moveTo(i);
    return #name#_stencilop(iter1_, iter2_); } 
									
  T_range_result operator()(const RectDomain<rank_>& d) const		
  { return T_range_result(iter1_(d), iter2_(d)); }			
									
  T_result operator[](int i) const					
  { return #name#_stencilop(iter1_[i], iter2_[i]); }				
									
  T_result fastRead(diffType i) const					
  {/* this probably isn't very fast... */				
    iter1_._bz_offsetData(i); iter2_._bz_offsetData(i);		
    T_result r = #name#_stencilop (iter1_, iter2_);				
    iter1_._bz_offsetData(-i); iter2_._bz_offsetData(-i);		
    return r;								
  }									

    /** This way of vectorizing won't work on stencils. */
    template<int N>
    typename tvresult<N>::Type fastRead_tv(diffType i) const {
      BZPRECHECK(0, "Can't vectorize stencils");
      return typename tvresult<N>::Type(iter1_.fastRead_tv<N>(i),
					iter2_.fastRead_tv<N>(i)); }
  
  T_result shift(int offset, int dim) const				
  {									
    iter1_._bz_offsetData(offset, dim);				
    iter2_._bz_offsetData(offset, dim);				
    T_result r = #name#_stencilop (iter1_, iter2_);				
    iter1_._bz_offsetData(-offset, dim);				
    iter2_._bz_offsetData(-offset, dim);				
    return r;								
  }									
									
  T_result shift(int offset1, int dim1, int offset2, int dim2) const	
  {									
    iter1_._bz_offsetData(offset1, dim1, offset2, dim2);		
    iter2_._bz_offsetData(offset1, dim1, offset2, dim2);		
    T_result r = #name#_stencilop (iter1_, iter2_);				
    iter1_._bz_offsetData(-offset1, dim1, -offset2, dim2);		
    iter2_._bz_offsetData(-offset1, dim1, -offset2, dim2);		
    return r;								
  }									
									
  void prettyPrint(BZ_STD_SCOPE(string) &str,				
		   prettyPrintFormat& format) const			
  {									
    str += "name (stencil)";						
    str += "(";							
    iter1_.prettyPrint(str, format);					
    str += ", ";							
    iter2_.prettyPrint(str, format);					
    str += ")";							
  }									
									
  template<typename T1, typename T2 = nilArraySection,		
	   class T3 = nilArraySection, typename T4 = nilArraySection,	
	   class T5 = nilArraySection, typename T6 = nilArraySection,	
	   class T7 = nilArraySection, typename T8 = nilArraySection,	
	   class T9 = nilArraySection, typename T10 = nilArraySection,	
	   class T11 = nilArraySection>					
  class SliceInfo {							
  public:								
    typedef #name#_et2<T_expr1, T_expr2, T_numtype> T_slice;	
  };								
    									
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
	   typename T7, typename T8, typename T9, typename T10, typename T11> 
  #name#_et2							
  operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const 
  {									
    return #name#_et2						
      (iter_(_bz_makeRange(r1),					
	     _bz_makeRange(r2),					
	     _bz_makeRange(r3),					
	     _bz_makeRange(r4),					
	     _bz_makeRange(r5),					
	     _bz_makeRange(r6),					
	     _bz_makeRange(r7),					
	     _bz_makeRange(r8),					
	     _bz_makeRange(r9),					
	     _bz_makeRange(r10),					
	     _bz_makeRange(r11)));					
  }									
};									
									
/* create ET object from application to expression */			
template<typename T1, typename T2>					
inline _bz_ArrayExpr<#name#_et2<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_range_result, #etresult#> > 
#name#(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1,				
     const BZ_BLITZ_SCOPE(ETBase)<T2>& d2)				
{									
  return _bz_ArrayExpr<#name#_et2<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_range_result, #etresult#> > 
    (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),MINB, MAXB)), 
     BZ_BLITZ_SCOPE(asExpr)<T2>::getExpr(d2.unwrap())(_bz_shrinkDomain(d2.unwrap().domain(),MINB, MAXB))); 
}									
/* matches to calls involving bare arrays (this is very annoying	
   because we have to exactly match every possible call combination	
   to ensure that this matches instead of the operator in		
   stencilops.h) */							
template<typename T1, typename T2, int N2>				
inline _bz_ArrayExpr<#name#_et2<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T2,N2> >::T_expr::T_range_result, #result#> > 
#name#(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1, Array<T2,N2>& d2)		
{ return #name#(d1.wrap(), d2.wrap()); }				
									
template<typename T1, typename T2, int N2>				
inline _bz_ArrayExpr<#name#_et2<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T2,N2> >::T_expr::T_range_result, #result#> > 
#name#(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1, const Array<T2,N2>& d2)	
{ return #name#(d1.wrap(), d2.wrap()); }				
									
template<typename T1, int N1, typename T2>				
inline _bz_ArrayExpr<#name#_et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T1,N1> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_range_result, #result#> > 
#name#(Array<T1,N1>& d1, const BZ_BLITZ_SCOPE(ETBase)<T2>& d2)		
{ return #name#(d1.wrap(), d2.wrap()); }				
									
template<typename T1, int N1, typename T2>				
inline _bz_ArrayExpr<#name#_et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T1,N1> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_range_result, #result#> > 
#name#(const Array<T1,N1>& d1, const BZ_BLITZ_SCOPE(ETBase)<T2>& d2)	
{ return #name#(d1.wrap(), d2.wrap()); }				
									
template<typename T1, int N1, typename T2, int N2>			
inline _bz_ArrayExpr<#name#_et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T1,N1> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T2,N2> >::T_expr::T_range_result, #result#> > 
#name#(const Array<T1,N1>& d1, Array<T2,N2>& d2)				
{ return #name#(d1.wrap(), d2.wrap()); }				
									
template<typename T1, int N1, typename T2, int N2>			
inline _bz_ArrayExpr<#name#_et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T1,N1> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T2,N2> >::T_expr::T_range_result, #result#> > 
#name#(Array<T1,N1>& d1, const Array<T2,N2>& d2)			
{ return #name#(d1.wrap(), d2.wrap()); }				
  									
template<typename T, int N>						
inline _bz_ArrayExpr<#name#_et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, #result#> > 
#name#(Array<T,N>& d1, Array<T,N>& d2)					
{ return #name#(d1.wrap(), d2.wrap()); }				
									
template<typename T, int N>						
inline _bz_ArrayExpr<#name#_et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, #result#> > 
#name#(const Array<T,N>& d1, const Array<T,N>& d2)			
{ return #name#(d1.wrap(), d2.wrap()); }"""
    os.write(stub.replace("#name#",name).replace("#result#",result).replace("#etresult#",etresult).replace("MINB",str(MINB)).replace("MAXB",str(MAXB)))

def BZ_ET_STENCILM(name, result_rank, MINB, MAXB):
    stub="""


/* Defines a stencil ET "#name#" that operates on an array
   array<P_numtype, N_rank> and returns a multicomponent
   array<TinyMatrix<P_numtype::T_element, rank, rank> >,
   N_rank>. P_numtype can be a TinyVector or a scalar, I think. */

template<typename P_expr>						
class #name#_et : public _bz_StencilExpr<P_expr, TinyMatrix<_bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element, result_rank, result_rank> > 
{									
public:									
  typedef _bz_StencilExpr<P_expr, TinyMatrix<_bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element, result_rank, result_rank> > T_base; 
  typedef _bz_typename T_base::T_numtype T_numtype;			
  typedef _bz_typename T_base::T_expr T_expr;				

  // there is no return type selection, as we are returning a
  // TinyMatrix. This must be returned as a FastTMCopyIterator since the
  // output of the stencil operator is a temporary.
  typedef ETBase<_bz_ArrayExpr<FastTM2CopyIterator<typename multicomponent_traits<typename P_expr::T_numtype>::T_element, result_rank, result_rank> > > T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;
  typedef T_numtype T_optype;

  template<int N> struct tvresult {
    typedef #name#_et<
      typename T_expr::template tvresult<N>::Type> Type;
  };

  typedef  #name#_et<_bz_typename P_expr::T_range_result> T_range_result; 
									 
  using T_base::iter_;						
  using T_base::rank_;						
public:								
  #name#_et(const #name#_et& a) :				
  _bz_StencilExpr<P_expr, T_numtype>(a)				
  { }								
     									
  #name#_et(BZ_ETPARM(T_expr) a) :					
  _bz_StencilExpr<P_expr, T_numtype>(a)				
  { }								
									 
  #name#_et(_bz_typename T_expr::T_ctorArg1 a) :			
  _bz_StencilExpr<P_expr, T_numtype>(a)				
  { }								
     									
  T_result operator*() const					
  { return #name#_stencilop(iter_); }						

  T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const 
  { iter_.moveTo(i); return #name#_stencilop(iter_); }				
     									
  T_range_result operator()(const RectDomain<rank_>& d) const		
  { return T_range_result(iter_(d)); }				
     									
  T_result operator[](int i) const					
  { return #name#_stencilop(iter_[i]); }						
     									
  T_result fastRead(diffType i) const				
  {/* this probably isn't very fast... */				
    iter_._bz_offsetData(i);						
    T_result r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-i);					
    return r;							
  }									

    /** This way of vectorizing won't work on stencils. */
    template<int N>
    typename tvresult<N>::Type fastRead_tv(diffType i) const {
      BZPRECHECK(0, "Can't vectorize stencils");
      return iter_.fastRead_tv<N>(i); }
     									
  T_result shift(int offset, int dim) const				
  {									
    iter_._bz_offsetData(offset, dim);				
    T_result r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset, dim);				
    return r;							
  }									
									 
  T_result shift(int offset1, int dim1, int offset2, int dim2) const 
  {									
    iter_._bz_offsetData(offset1, dim1, offset2, dim2);		
    T_result r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		
    return r;							
  }									
									 
  void prettyPrint(BZ_STD_SCOPE(string) &str,			
		   prettyPrintFormat& format) const			
  {									
    str += "name (stencil)";						
    str += "(";							
    iter_.prettyPrint(str, format);					
    str += ")";							
  }									
									
  template<typename T1, typename T2 = nilArraySection,		
	   class T3 = nilArraySection, typename T4 = nilArraySection,	
	   class T5 = nilArraySection, typename T6 = nilArraySection,	
	   class T7 = nilArraySection, typename T8 = nilArraySection,	
	   class T9 = nilArraySection, typename T10 = nilArraySection,	
	   class T11 = nilArraySection>					
  class SliceInfo {							
  public:								
    typedef #name#_et<T_expr> T_slice;				
  };								
    									
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
	   typename T7, typename T8, typename T9, typename T10, typename T11> 
  #name#_et							
  operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const 
  {									
    return #name#_et						
      (iter_(_bz_makeRange(r1),					
	     _bz_makeRange(r2),					
	     _bz_makeRange(r3),					
	     _bz_makeRange(r4),					
	     _bz_makeRange(r5),					
	     _bz_makeRange(r6),					
	     _bz_makeRange(r7),					
	     _bz_makeRange(r8),					
	     _bz_makeRange(r9),					
	     _bz_makeRange(r10),					
	     _bz_makeRange(r11)));					
  }									
};									
/* create ET from application to expression */			
template<typename T1>						
inline _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> > 
#name#(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1)				
{									
  return _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	
    (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),MINB, MAXB))); 
}								
"""
    os.write(stub.replace("#name#",name).replace("result_rank",str(result_rank)).replace("MINB",str(MINB)).replace("MAXB",str(MAXB)))
    BZ_ET_STENCIL_REDIRECT(name)						

def BZ_ET_STENCILV(name, result_rank, MINB, MAXB):
    stub="""


/* Defines a stencil ET "#name#" that operates on a (scalar) array<P_numtype,
   N_rank> and returns a multicomponent
   array<TinyVector<P_numtype::T_element, result_rank> >, N_rank>. */

template<typename P_expr>						
class #name#_et : public _bz_StencilExpr<P_expr, TinyVector<typename P_expr::T_numtype,result_rank> > 
{									
public:									
  typedef _bz_StencilExpr<P_expr, TinyVector<typename P_expr::T_numtype,result_rank> > T_base; 
  typedef _bz_typename T_base::T_numtype T_numtype;			
  typedef _bz_typename T_base::T_expr T_expr;				

  // there is no return type selection, we assume P_numtype is scalar
  // and that we are returning a TinyVector. This needs to be returned
  // as a FastTVCopyIterator that keeps a copy of the TV it is
  // iterating over, since the result of the stencil operator is a temporary.
  typedef ETBase<_bz_ArrayExpr<FastTV2CopyIterator<typename P_expr::T_numtype, result_rank> > > T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;
  typedef typename T_expr::T_numtype T_optype;

  template<int N> struct tvresult {
    typedef #name#_et<
      typename T_expr::template tvresult<N>::Type> Type;
  };

  typedef  #name#_et<_bz_typename P_expr::T_range_result> T_range_result; 
									 
  using T_base::iter_;						
  using T_base::rank_;							
public:								
  #name#_et(const #name#_et& a) :					
  _bz_StencilExpr<P_expr, T_numtype>(a)				
  { }								
									 
  #name#_et(BZ_ETPARM(T_expr) a) :					
  _bz_StencilExpr<P_expr, T_numtype>(a)				
  { }								
									 
  #name#_et(_bz_typename T_expr::T_ctorArg1 a) :			
  _bz_StencilExpr<P_expr, T_numtype>(a)				
  { }								
									 
  T_result operator*() const					
  { return #name#_stencilop(iter_); }						

  T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const 
  { iter_.moveTo(i); return #name#_stencilop(iter_); }				
									 
  T_range_result operator()(const RectDomain<rank_>& d) const		
  { return T_range_result(iter_(d)); }				
									 
  T_result operator[](int i) const					
  { return #name#_stencilop(iter_[i]); }						
									 
  T_result fastRead(diffType i) const				
  {/* this probably isn't very fast... */				
    iter_._bz_offsetData(i);						
    T_numtype r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-i);					
    return r;							
  }									
									 
    /** This way of vectorizing won't work on stencils. */
    template<int N>
    typename tvresult<N>::Type fastRead_tv(diffType i) const {
      BZPRECHECK(0, "Can't vectorize stencils");
      return iter_.fastRead_tv<N>(i); }

  T_numtype shift(int offset, int dim) const				
  {									
    iter_._bz_offsetData(offset, dim);				
    T_numtype r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset, dim);				
    return r;								
  }									
									 
  T_numtype shift(int offset1, int dim1, int offset2, int dim2) const 
  {									
    iter_._bz_offsetData(offset1, dim1, offset2, dim2);		
    T_numtype r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		
    return r;								
  }									
									 
  void prettyPrint(BZ_STD_SCOPE(string) &str,				
		   prettyPrintFormat& format) const			
  {									
    str += "name (stencil)";						
    str += "(";							
    iter_.prettyPrint(str, format);					
    str += ")";							
  }									
									
  template<typename T1, typename T2 = nilArraySection,		
	   class T3 = nilArraySection, typename T4 = nilArraySection,	
	   class T5 = nilArraySection, typename T6 = nilArraySection,	
	   class T7 = nilArraySection, typename T8 = nilArraySection,	
	   class T9 = nilArraySection, typename T10 = nilArraySection,	
	   class T11 = nilArraySection>					
  class SliceInfo {							
  public:								
    typedef #name#_et<T_expr> T_slice;				
  };								
    									
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
	   typename T7, typename T8, typename T9, typename T10, typename T11> 
  #name#_et							
  operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const 
  {									
    return #name#_et						
      (iter_(_bz_makeRange(r1),					
	     _bz_makeRange(r2),					
	     _bz_makeRange(r3),					
	     _bz_makeRange(r4),					
	     _bz_makeRange(r5),					
	     _bz_makeRange(r6),					
	     _bz_makeRange(r7),					
	     _bz_makeRange(r8),					
	     _bz_makeRange(r9),					
	     _bz_makeRange(r10),					
	     _bz_makeRange(r11)));					
  }									
};									
/* create ET from application to expression */				
template<typename T1>							
inline _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	
#name#(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1)				
{									
  return _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	
    (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),MINB, MAXB))); 
}									
"""
    os.write(stub.replace("#name#",name).replace("result_rank",str(result_rank)).replace("MINB",str(MINB)).replace("MAXB",str(MAXB)))
    BZ_ET_STENCIL_REDIRECT(name)						


def BZ_ET_STENCIL_SCA(name, MINB, MAXB):
    stub="""


/** Defines a stencil ET "#name#" that operates on a multicomponent
   array<P_numtype, N_rank> and returns a scalar
   array<P_numtype::T_element, N_rank>. */

template<typename P_expr>						
class #name#_et : public _bz_StencilExpr<P_expr, _bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element> 
{									
public:									
  typedef _bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element T_result; 
  typedef _bz_StencilExpr<P_expr, T_result> T_base;			
  typedef _bz_typename T_base::T_numtype T_numtype;			
  typedef _bz_typename T_base::T_expr T_expr;				

  // there is no selecting return type here. because we *know* it is
  // scalar T_result, there's no question of whether we could be doing
  // multicomponent evaluations.    
  typedef T_result T_typeprop;
  typedef T_numtype T_optype;

  template<int N> struct tvresult {
    typedef #name#_et<
      typename T_expr::template tvresult<N>::Type> Type;
  };

typedef  #name#_et<_bz_typename P_expr::T_range_result> T_range_result; 
									 
  using T_base::iter_;						
  using T_base::rank_;							
public:								
  #name#_et(const #name#_et& a) :					
  _bz_StencilExpr<P_expr, T_numtype>(a)				
  { }								
									 
  #name#_et(BZ_ETPARM(T_expr) a) :					
  _bz_StencilExpr<P_expr, T_numtype>(a)				
  { }								
									 
  #name#_et(_bz_typename T_expr::T_ctorArg1 a) :			
  _bz_StencilExpr<P_expr, T_numtype>(a)				
  { }								
    									
  T_numtype operator*() const						
  { return #name#_stencilop(iter_); }						
  T_numtype operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const 
  { iter_.moveTo(i); return #name#_stencilop(iter_); }				
									
  T_range_result operator()(const RectDomain<rank_>& d) const		
  { return T_range_result(iter_(d)); }				
									
  T_numtype operator[](int i) const					
  { return #name#_stencilop(iter_[i]); }						
									
  T_numtype fastRead(diffType i) const				
  {/* this probably isn't very fast... */				
    iter_._bz_offsetData(i);						
    T_numtype r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-i);						
    return r;								
  }									
									
    /** This way of vectorizing won't work on stencils. */
    template<int N>
    typename tvresult<N>::Type fastRead_tv(diffType i) const {
      BZPRECHECK(0, "Can't vectorize stencils");
      return iter_.fastRead_tv<N>(i); }

  T_numtype shift(int offset, int dim) const				
  {									
    iter_._bz_offsetData(offset, dim);				
    T_numtype r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset, dim);				
    return r;								
  }									
									
  T_numtype shift(int offset1, int dim1, int offset2, int dim2) const	
  {									
    iter_._bz_offsetData(offset1, dim1, offset2, dim2);		
    T_numtype r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		
    return r;								
  }									
									
  void prettyPrint(BZ_STD_SCOPE(string) &str,				
		   prettyPrintFormat& format) const			
  {									
    str += "name (stencil)";						
    str += "(";							
    iter_.prettyPrint(str, format);					
    str += ")";							
  }									
									
  template<typename T1, typename T2 = nilArraySection,		
	   class T3 = nilArraySection, typename T4 = nilArraySection,	
	   class T5 = nilArraySection, typename T6 = nilArraySection,	
	   class T7 = nilArraySection, typename T8 = nilArraySection,	
	   class T9 = nilArraySection, typename T10 = nilArraySection,	
	   class T11 = nilArraySection>					
  class SliceInfo {							
  public:								
    typedef #name#_et<T_expr> T_slice;				
  };								
    									
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
	   typename T7, typename T8, typename T9, typename T10, typename T11> 
  #name#_et							
  operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const 
  {									
    return #name#_et						
      (iter_(_bz_makeRange(r1),					
	     _bz_makeRange(r2),					
	     _bz_makeRange(r3),					
	     _bz_makeRange(r4),					
	     _bz_makeRange(r5),					
	     _bz_makeRange(r6),					
	     _bz_makeRange(r7),					
	     _bz_makeRange(r8),					
	     _bz_makeRange(r9),					
	     _bz_makeRange(r10),					
	     _bz_makeRange(r11)));					
  }									
};									
/* create ET from application to expression */			
template<typename T1>							
inline _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	
#name#(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1)				
{									
  return _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	
    (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),MINB, MAXB))); 
}									
"""
    os.write(stub.replace("#name#",name).replace("MINB",str(MINB)).replace("MAXB",str(MAXB)))
    BZ_ET_STENCIL_REDIRECT(name)						


def BZ_ET_STENCIL_DIFF(name, MINB, MAXB):
    stub="""


/* Defines a stencil ET difference operator "#name#" that operates on
   an array<P_numtype, N_rank> and returns an array of identical
   type. (The only significance of the "difference" aspect is that the
   operator is assumed to take a second argument which is the
   dimension to do the difference in). MINB and MAXB are integer
   expressions describing the extent of the operator in the operating
   dimension. */

template<typename P_expr>						
class #name#_et : public _bz_StencilExpr<P_expr, _bz_typename P_expr::T_numtype>	
{									
public:								
  typedef _bz_StencilExpr<P_expr, _bz_typename P_expr::T_numtype> T_base; 
  typedef _bz_typename T_base::T_numtype T_numtype;			
  typedef _bz_typename T_base::T_expr T_expr;				

    // select return type
  typedef typename unwrapET<typename T_expr::T_result>::T_unwrapped test;
  typedef typename selectET<typename T_expr::T_typeprop, 
			    T_numtype, 
			    #name#_et<test> >::T_selected T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;
  typedef T_numtype T_optype;

  template<int N> struct tvresult {
    typedef #name#_et<
      typename T_expr::template tvresult<N>::Type> Type;
  };

typedef  #name#_et<_bz_typename P_expr::T_range_result> T_range_result; 
									
  using T_base::iter_;						
  using T_base::rank_;							
public:								
  #name#_et(const #name#_et& a) :				
  _bz_StencilExpr<P_expr, T_numtype>(a), dim_(a.dim_)	
  { }								
									
  #name#_et(BZ_ETPARM(T_expr) a, int dim) :			
  _bz_StencilExpr<P_expr, T_numtype>(a), dim_(dim)		
  { }								
									
  #name#_et(_bz_typename T_expr::T_ctorArg1 a, int dim) :	
  _bz_StencilExpr<P_expr, T_numtype>(a), dim_(dim)		
  { }								
									
  T_result operator*() const						
  { return #name#_stencilop(iter_, dim_); }					
  T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const 
  { iter_.moveTo(i); return #name#_stencilop(iter_, dim_); }			
									
  T_range_result operator()(const RectDomain<rank_>& d) const		
  { return T_range_result(iter_(d), dim_); }				
									
  T_result operator[](int i) const					
  { return #name#_stencilop(iter_[i], dim_); }					
									
  T_result fastRead(diffType i) const				
  {/* this probably isn't very fast... */				
    iter_._bz_offsetData(i);						
    T_result r = #name#_stencilop (iter_, dim_);					
    iter_._bz_offsetData(-i);						
    return r;								
  }									
									
    /** This way of vectorizing won't work on stencils. */
    template<int N>
    typename tvresult<N>::Type fastRead_tv(diffType i) const {
      BZPRECHECK(0, "Can't vectorize stencils");
      return typename tvresult<N>::Type(iter_.fastRead_tv<N>(i),dim_); }

  T_result shift(int offset, int dim) const				
  {									
    iter_._bz_offsetData(offset, dim);				
    T_result r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset, dim);				
    return r;								
  }									
									
  T_result shift(int offset1, int dim1, int offset2, int dim2) const	
  {									
    iter_._bz_offsetData(offset1, dim1, offset2, dim2);		
    T_result r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		
    return r;								
  }									
									
  void prettyPrint(BZ_STD_SCOPE(string) &str,				
		   prettyPrintFormat& format) const			
  {									
    str += "name (stencil)";						
    str += "(";							
    iter_.prettyPrint(str, format);					
    str += ")";							
  }									
									
  template<typename T1, typename T2 = nilArraySection,		
	   class T3 = nilArraySection, typename T4 = nilArraySection,	
	   class T5 = nilArraySection, typename T6 = nilArraySection,	
	   class T7 = nilArraySection, typename T8 = nilArraySection,	
	   class T9 = nilArraySection, typename T10 = nilArraySection,	
	   class T11 = nilArraySection>					
  class SliceInfo {							
  public:								
    typedef #name#_et<T_expr> T_slice;				
  };								
    									
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
	   typename T7, typename T8, typename T9, typename T10, typename T11> 
  #name#_et							
  operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const 
  {									
    return #name#_et						
      (iter_(_bz_makeRange(r1),					
	     _bz_makeRange(r2),					
	     _bz_makeRange(r3),					
	     _bz_makeRange(r4),					
	     _bz_makeRange(r5),					
	     _bz_makeRange(r6),					
	     _bz_makeRange(r7),					
	     _bz_makeRange(r8),					
	     _bz_makeRange(r9),					
	     _bz_makeRange(r10),					
	     _bz_makeRange(r11)),dim_);				
  }									
									
private:								
  int dim_;								
};									
/* create ET from application to expression */				
template<typename T1>							
inline _bz_ArrayExpr<#name#_et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	
#name#(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1, int dim)			
{									
  TinyVector<int, BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::rank_> minb(0), maxb(0); 
  minb[dim]=MINB; maxb[dim]=MAXB;					
  return _bz_ArrayExpr<#name#_et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	
    (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),minb, maxb)), dim); 
}									
/* forward operations on arrays to main function */			
template<typename T, int N>						
inline _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > 
#name#(const Array<T,N>& d1, int dim)					
{ return #name#(d1.wrap(), dim); }					
template<typename T, int N>						
inline _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > 
#name#(Array<T,N>& d1, int dim)						
{ return #name#(d1.wrap(), dim); }
"""
    os.write(stub.replace("#name#",name).replace("MINB",str(MINB)).replace("MAXB",str(MAXB)))


def BZ_ET_STENCIL_MULTIDIFF(name, MINB, MAXB):
    stub="""


/** Defines a stencil ET difference operator "#name#" that operates on a
   multicomponent array<P_numtype, N_rank> and returns an
   array<P_numtype::T_element, N_rank>. */
template<typename P_expr>						
class #name#_et_multi : public _bz_StencilExpr<P_expr, _bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element> 
{									
public:								
  typedef _bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element T_result; 
  typedef _bz_StencilExpr<P_expr, T_result> T_base;		
  typedef _bz_typename T_base::T_numtype T_numtype;			
  typedef _bz_typename T_base::T_expr T_expr;				

  // there is no selecting return type here. because we *know* it is
  // T_result, there's no question of whether we could be doing
   // multicomponent evaluations.    
  typedef T_result T_typeprop;
  typedef T_numtype T_optype;

  template<int N> struct tvresult {
    typedef #name#_et_multi<
      typename T_expr::template tvresult<N>::Type> Type;
  };

typedef  #name#_et_multi<_bz_typename P_expr::T_range_result> T_range_result; 
									
  using T_base::iter_;						
  using T_base::rank_;							
public:								
  #name#_et_multi(const #name#_et_multi& a) :		
  _bz_StencilExpr<P_expr, T_numtype>(a), comp_(a.comp_), dim_(a.dim_) 
  { }								

  #name#_et_multi(BZ_ETPARM(T_expr) a, int comp, int dim) :	
  _bz_StencilExpr<P_expr, T_numtype>(a),			
    comp_(comp), dim_(dim)						
  { }								
									
  #name#_et_multi(_bz_typename T_expr::T_ctorArg1 a, int comp, int dim) : 
  _bz_StencilExpr<P_expr, T_numtype>(a),			
    comp_(comp), dim_(dim)						
  { }								
									
  T_numtype operator*() const						
  { return #name#_stencilop(iter_, comp_, dim_); }				
  T_numtype operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const 
  { iter_.moveTo(i); return #name#_stencilop(iter_, comp_, dim_); }		
									
  T_range_result operator()(const RectDomain<rank_>& d) const		
  { return T_range_result(iter_(d), comp_, dim_); }			
									
  T_numtype operator[](int i) const					
  { return #name#_stencilop(iter_[i], comp_, dim_); }				
									
  T_numtype fastRead(diffType i) const				
  {/* this probably isn't very fast... */				
    iter_._bz_offsetData(i);						
    T_numtype r = #name#_stencilop (iter_, comp_, dim_);				
    iter_._bz_offsetData(-i);						
    return r;								
  }									
									
    /** This way of vectorizing won't work on stencils. */
    template<int N>
    typename tvresult<N>::Type fastRead_tv(diffType i) const {
      BZPRECHECK(0, "Can't vectorize stencils");
      return typename tvresult<N>::Type(iter_.fastRead_tv<N>(i),comp_,dim_); }

  T_numtype shift(int offset, int dim) const				
  {									
    iter_._bz_offsetData(offset, dim);				
    T_numtype r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset, dim);				
    return r;								
  }									
									
  T_numtype shift(int offset1, int dim1, int offset2, int dim2) const	
  {									
    iter_._bz_offsetData(offset1, dim1, offset2, dim2);		
    T_numtype r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		
    return r;								
  }									
									
  void prettyPrint(BZ_STD_SCOPE(string) &str,				
		   prettyPrintFormat& format) const			
  {									
    str += "name (stencil)";						
    str += "(";							
    iter_.prettyPrint(str, format);					
    str += ")";							
  }									
									
  template<typename T1, typename T2 = nilArraySection,		
	   class T3 = nilArraySection, typename T4 = nilArraySection,	
	   class T5 = nilArraySection, typename T6 = nilArraySection,	
	   class T7 = nilArraySection, typename T8 = nilArraySection,	
	   class T9 = nilArraySection, typename T10 = nilArraySection,	
	   class T11 = nilArraySection>					
  class SliceInfo {							
  public:								
    typedef #name#_et_multi<T_expr> T_slice;				
  };								
    									
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
	   typename T7, typename T8, typename T9, typename T10, typename T11> 
  #name#_et_multi							
  operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const 
  {									
    return #name#_et_multi						
      (iter_(_bz_makeRange(r1),					
	     _bz_makeRange(r2),					
	     _bz_makeRange(r3),					
	     _bz_makeRange(r4),					
	     _bz_makeRange(r5),					
	     _bz_makeRange(r6),					
	     _bz_makeRange(r7),					
	     _bz_makeRange(r8),					
	     _bz_makeRange(r9),					
	     _bz_makeRange(r10),					
	     _bz_makeRange(r11)),comp_, dim_);			
  }									
									
private:								
  int comp_;								
  int dim_;								
};									
/* create ET from application to expression */			
template<typename T1>							
inline _bz_ArrayExpr<#name#_et_multi<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	
#name#(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1, int comp, int dim)		
{									
  TinyVector<int, BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::rank_> minb(0), maxb(0);	
  minb[dim]=MINB; maxb[dim]=MAXB;					
  return _bz_ArrayExpr<#name#_et_multi<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> > 
    (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),minb, maxb)), comp, dim); 
}									
/* forward operations on arrays to main function */			
template<typename T, int N>						
inline _bz_ArrayExpr<#name#_et_multi<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > 
#name#(const Array<T,N>& d1, int comp, int dim)				
{ return #name#(d1.wrap(), comp, dim); }				
									
template<typename T, int N>						
inline _bz_ArrayExpr<#name#_et_multi<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > 
#name#(Array<T,N>& d1, int comp, int dim)				
{ return #name#(d1.wrap(), comp, dim); }
"""
    os.write(stub.replace("#name#",name).replace("MINB",str(MINB)).replace("MAXB",str(MAXB)))
    BZ_ET_STENCIL_REDIRECT(name)						


def BZ_ET_STENCIL_DIFF2(name, MINB1, MAXB1, MINB2, MAXB2):
    stub="""


/** Defines a stencil ET double-difference operator "#name#" that
   operates on an array<P_numtype, N_rank> and returns an array of
   identical type. (The only significance of the "double-difference" aspect
   is that the operator is assumed to take two extra arguments which
   are the dimensions to do the differences in). */

template<typename P_expr>						
class #name#_et : public _bz_StencilExpr<P_expr, _bz_typename P_expr::T_numtype> 
{									
public:								
  typedef _bz_StencilExpr<P_expr, _bz_typename P_expr::T_numtype> T_base;	
  typedef _bz_typename T_base::T_numtype T_numtype;			
  typedef _bz_typename T_base::T_expr T_expr;				

    // select return type
  typedef typename unwrapET<typename T_expr::T_result>::T_unwrapped test;
  typedef typename selectET<typename T_expr::T_typeprop, 
			    T_numtype, 
			    #name#_et<test> >::T_selected T_typeprop;
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;
  typedef T_numtype T_optype;

  template<int N> struct tvresult {
    typedef #name#_et<
      typename T_expr::template tvresult<N>::Type> Type;
  };

typedef  #name#_et<_bz_typename P_expr::T_range_result> T_range_result; 
   									
  using T_base::iter_;							
  using T_base::rank_;							
public:								
  #name#_et(const #name#_et& a) :					
  _bz_StencilExpr<P_expr, T_numtype>(a),				
    dim1_(a.dim1_), dim2_(a.dim2_)					
  { }								
   									
  #name#_et(BZ_ETPARM(T_expr) a, int dim1, int dim2) :		
  _bz_StencilExpr<P_expr, T_numtype>(a),				
    dim1_(dim1), dim2_(dim2)						
  { }									
   									
  #name#_et(_bz_typename T_expr::T_ctorArg1 a,			
	      int dim1, int dim2) :					
  _bz_StencilExpr<P_expr, T_numtype>(a),				
    dim1_(dim1), dim2_(dim2)						
  { }									
   									
  T_numtype operator*() const						
  { return #name#_stencilop(iter_, dim1_, dim2_); }				
  T_numtype operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const 
  { iter_.moveTo(i); return #name#_stencilop(iter_, dim1_, dim2_); }		
									
  T_range_result operator()(const RectDomain<rank_>& d) const		
  { return T_range_result(iter_(d), dim1_, dim2_); }			
   									
  T_numtype operator[](int i) const					
  { return #name#_stencilop(iter_[i], dim1_, dim2_); }				
									
  T_numtype fastRead(diffType i) const					
  {/* this probably isn't very fast... */				
    iter_._bz_offsetData(i);						
    T_numtype r = #name#_stencilop (iter_, dim1_, dim2_);				
    iter_._bz_offsetData(-i);						
    return r;								
  }									
									
    /** This way of vectorizing won't work on stencils. */
    template<int N>
    typename tvresult<N>::Type fastRead_tv(diffType i) const {
      BZPRECHECK(0, "Can't vectorize stencils");
      return typename tvresult<N>::Type(iter_.fastRead_tv<N>(i),dim1_,dim2_); }

  T_numtype shift(int offset, int dim) const				
  {									
    iter_._bz_offsetData(offset, dim);				
    T_numtype r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset, dim);				
    return r;								
  }									
									
  T_numtype shift(int offset1, int dim1, int offset2, int dim2) const	
  {									
    iter_._bz_offsetData(offset1, dim1, offset2, dim2);		
    T_numtype r = #name#_stencilop (iter_);					
    iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		
    return r;								
  }									
									
  void prettyPrint(BZ_STD_SCOPE(string) &str,				
		   prettyPrintFormat& format) const			
  {									
    str += "name (stencil)";						
    str += "(";							
    iter_.prettyPrint(str, format);					
    str += ")";							
  }									
									
  template<typename T1, typename T2 = nilArraySection,		
	   class T3 = nilArraySection, typename T4 = nilArraySection,	
	   class T5 = nilArraySection, typename T6 = nilArraySection,	
	   class T7 = nilArraySection, typename T8 = nilArraySection,	
	   class T9 = nilArraySection, typename T10 = nilArraySection,	
	   class T11 = nilArraySection>					
  class SliceInfo {							
  public:								
    typedef #name#_et<T_expr> T_slice;				
  };								
    									
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, 
	   typename T7, typename T8, typename T9, typename T10, typename T11> 
  #name#_et							
  operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const 
  {									
    return #name#_et						
      (iter_(_bz_makeRange(r1),					
	     _bz_makeRange(r2),					
	     _bz_makeRange(r3),					
	     _bz_makeRange(r4),					
	     _bz_makeRange(r5),					
	     _bz_makeRange(r6),					
	     _bz_makeRange(r7),					
	     _bz_makeRange(r8),					
	     _bz_makeRange(r9),					
	     _bz_makeRange(r10),					
	     _bz_makeRange(r11)), dim1_, dim2_);			
  }									
   									
private:								
  int dim1_, dim2_;							
};									
 									
/* create ET from application to expression */
template<typename T1>							
inline _bz_ArrayExpr<#name#_et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	
#name#(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1, int dim1, int dim2)		
{									
  TinyVector<int, BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::rank_> minb(0), maxb(0);	
  minb[dim1]=MINB1; maxb[dim1]=MAXB1;					
  minb[dim2]=MINB2; maxb[dim2]=MAXB2;					
  return _bz_ArrayExpr<#name#_et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	
    (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),minb, maxb)), dim1, dim2); 
}									
/* forward operations on arrays to main function */			
template<typename T, int N>						
inline _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > 
#name#(const Array<T,N>& d1, int dim1, int dim2)			
{ return #name#(d1.wrap(), dim1, dim2); }				
									
template<typename T, int N>						
inline _bz_ArrayExpr<#name#_et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > 
#name#(Array<T,N>& d1, int dim1, int dim2)				
{ return #name#(d1.wrap(), dim1, dim2); }
"""
    os.write(stub.replace("#name#",name).replace("MINB1",str(MINB1)).replace("MAXB1",str(MAXB1)).replace("MINB2",str(MINB2)).replace("MAXB2",str(MAXB2)))


# now generate the file

os.write("""
/* Definitions of ET stencil operators. This file was generated by
   generate/genstencils.py. Do not edit.

   Note: You can't pass templates with >1 parameter as macro
   parameters because cpp doesn't recognize that the comma is balanced
   between the angle brackets and interprets them as multiple
   arguments, i.e., the following alternative declaration of grad2D
   will not work:

   BZ_ET_STENCIL(grad2D, TinyVector<P_numtype, 2>, 
   TinyVector<typename T1::T_numtype, 2>, shape(-1,1), shape(1,1))

   instead, you have to use the above bzCC ("ConCatenate") macro to
   protect the things containing commas. The following would work:

   BZ_ET_STENCIL(grad2D, bzCC(TinyVector<P_numtype, 2>), 
   bzCC(TinyVector<typename T1::T_numtype, 2>), shape(-1,-1), shape(1,1))

*/

BZ_NAMESPACE(blitz)

""")

BZ_ET_STENCIL_DIFF("central12", -1, 1)
BZ_ET_STENCIL_DIFF("central22", -1, 1)
BZ_ET_STENCIL_DIFF("central32", -2, 2)
BZ_ET_STENCIL_DIFF("central42", -2, 2)
BZ_ET_STENCIL_DIFF("central14", -2, 2)
BZ_ET_STENCIL_DIFF("central24", -2, 2)
BZ_ET_STENCIL_DIFF("central34", -2, 2)
BZ_ET_STENCIL_DIFF("central44", -2, 2)
BZ_ET_STENCIL_DIFF("central12n", -1, 1)
BZ_ET_STENCIL_DIFF("central22n", -1, 1)
BZ_ET_STENCIL_DIFF("central32n", -2, 2)
BZ_ET_STENCIL_DIFF("central42n", -2, 2)
BZ_ET_STENCIL_DIFF("central14n", -2, 2)
BZ_ET_STENCIL_DIFF("central24n", -2, 2)
BZ_ET_STENCIL_DIFF("central34n", -2, 2)
BZ_ET_STENCIL_DIFF("central44n", -2, 2)

BZ_ET_STENCIL_DIFF("backward11", -1, 0)
BZ_ET_STENCIL_DIFF("backward21", -2, 0)
BZ_ET_STENCIL_DIFF("backward31", -3, 0)
BZ_ET_STENCIL_DIFF("backward41", -4, 0)
BZ_ET_STENCIL_DIFF("backward12", -2, 0)
BZ_ET_STENCIL_DIFF("backward22", -3, 0)
BZ_ET_STENCIL_DIFF("backward32", -4, 0)
BZ_ET_STENCIL_DIFF("backward42", -5, 0)
BZ_ET_STENCIL_DIFF("backward11n", -1, 0)
BZ_ET_STENCIL_DIFF("backward21n", -2, 0)
BZ_ET_STENCIL_DIFF("backward31n", -3, 0)
BZ_ET_STENCIL_DIFF("backward41n", -4, 0)
BZ_ET_STENCIL_DIFF("backward12n", -2, 0)
BZ_ET_STENCIL_DIFF("backward22n", -3, 0)
BZ_ET_STENCIL_DIFF("backward32n", -4, 0)
BZ_ET_STENCIL_DIFF("backward42n", -5, 0)

BZ_ET_STENCIL_DIFF("forward11", 0, 1)
BZ_ET_STENCIL_DIFF("forward21", 0, 2)
BZ_ET_STENCIL_DIFF("forward31", 0, 3)
BZ_ET_STENCIL_DIFF("forward41", 0, 4)
BZ_ET_STENCIL_DIFF("forward12", 0, 2)
BZ_ET_STENCIL_DIFF("forward22", 0, 3)
BZ_ET_STENCIL_DIFF("forward32", 0, 4)
BZ_ET_STENCIL_DIFF("forward42", 0, 5)
BZ_ET_STENCIL_DIFF("forward11n", 0, 1)
BZ_ET_STENCIL_DIFF("forward21n", 0, 2)
BZ_ET_STENCIL_DIFF("forward31n", 0, 3)
BZ_ET_STENCIL_DIFF("forward41n", 0, 4)
BZ_ET_STENCIL_DIFF("forward12n", 0, 2)
BZ_ET_STENCIL_DIFF("forward22n", 0, 3)
BZ_ET_STENCIL_DIFF("forward32n", 0, 4)
BZ_ET_STENCIL_DIFF("forward42n", 0, 5)

BZ_ET_STENCIL_MULTIDIFF("central12", -1, 1)
BZ_ET_STENCIL_MULTIDIFF("central22", -1, 1)
BZ_ET_STENCIL_MULTIDIFF("central32", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central42", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central14", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central24", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central34", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central44", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central12n", -1, 1)
BZ_ET_STENCIL_MULTIDIFF("central22n", -1, 1)
BZ_ET_STENCIL_MULTIDIFF("central32n", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central42n", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central14n", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central24n", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central34n", -2, 2)
BZ_ET_STENCIL_MULTIDIFF("central44n", -2, 2)

BZ_ET_STENCIL_MULTIDIFF("backward11", -1, 0)
BZ_ET_STENCIL_MULTIDIFF("backward21", -2, 0)
BZ_ET_STENCIL_MULTIDIFF("backward31", -3, 0)
BZ_ET_STENCIL_MULTIDIFF("backward41", -4, 0)
BZ_ET_STENCIL_MULTIDIFF("backward12", -2, 0)
BZ_ET_STENCIL_MULTIDIFF("backward22", -3, 0)
BZ_ET_STENCIL_MULTIDIFF("backward32", -4, 0)
BZ_ET_STENCIL_MULTIDIFF("backward42", -5, 0)
BZ_ET_STENCIL_MULTIDIFF("backward11n", -1, 0)
BZ_ET_STENCIL_MULTIDIFF("backward21n", -2, 0)
BZ_ET_STENCIL_MULTIDIFF("backward31n", -3, 0)
BZ_ET_STENCIL_MULTIDIFF("backward41n", -4, 0)
BZ_ET_STENCIL_MULTIDIFF("backward12n", -2, 0)
BZ_ET_STENCIL_MULTIDIFF("backward22n", -3, 0)
BZ_ET_STENCIL_MULTIDIFF("backward32n", -4, 0)
BZ_ET_STENCIL_MULTIDIFF("backward42n", -5, 0)

BZ_ET_STENCIL_MULTIDIFF("forward11", 0, 1)
BZ_ET_STENCIL_MULTIDIFF("forward21", 0, 2)
BZ_ET_STENCIL_MULTIDIFF("forward31", 0, 3)
BZ_ET_STENCIL_MULTIDIFF("forward41", 0, 4)
BZ_ET_STENCIL_MULTIDIFF("forward12", 0, 2)
BZ_ET_STENCIL_MULTIDIFF("forward22", 0, 3)
BZ_ET_STENCIL_MULTIDIFF("forward32", 0, 4)
BZ_ET_STENCIL_MULTIDIFF("forward42", 0, 5)
BZ_ET_STENCIL_MULTIDIFF("forward11n", 0, 1)
BZ_ET_STENCIL_MULTIDIFF("forward21n", 0, 2)
BZ_ET_STENCIL_MULTIDIFF("forward31n", 0, 3)
BZ_ET_STENCIL_MULTIDIFF("forward41n", 0, 4)
BZ_ET_STENCIL_MULTIDIFF("forward12n", 0, 2)
BZ_ET_STENCIL_MULTIDIFF("forward22n", 0, 3)
BZ_ET_STENCIL_MULTIDIFF("forward32n", 0, 4)
BZ_ET_STENCIL_MULTIDIFF("forward42n", 0, 5)

BZ_ET_STENCIL("Laplacian2D", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCIL("Laplacian3D", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCIL("Laplacian2D4", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-2,-2)", "shape(2,2)")
BZ_ET_STENCIL("Laplacian2D4n", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-2,-2)", "shape(2,2)")
BZ_ET_STENCIL("Laplacian3D4", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-2,-2,-2)", "shape(2,2,2)")
BZ_ET_STENCIL("Laplacian3D4n", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-2,-2,-2)", "shape(2,2,2)")

BZ_ET_STENCILV("grad2D", 2, "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCILV("grad2D4", 2, "shape(-2,-2)", "shape(2,2)")
BZ_ET_STENCILV("grad3D", 3, "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCILV("grad3D4", 3, "shape(-2,-2,-2)", "shape(2,2,2)")
BZ_ET_STENCILV("grad2Dn", 2, "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCILV("grad2D4n", 2, "shape(-2,-2)", "shape(2,2)")
BZ_ET_STENCILV("grad3Dn", 3, "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCILV("grad3D4n", 3, "shape(-2,-2,-2)", "shape(2,2,2)")
BZ_ET_STENCILV("gradSqr2D", 2, "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCILV("gradSqr2D4", 2, "shape(-2,-2)", "shape(2,2)")
BZ_ET_STENCILV("gradSqr3D", 3, "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCILV("gradSqr3D4", 3, "shape(-2,-2,-2)", "shape(2,2,2)")
BZ_ET_STENCILV("gradSqr2Dn", 2, "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCILV("gradSqr2D4n", 2, "shape(-2,-2)", "shape(2,2)")
BZ_ET_STENCILV("gradSqr3Dn", 3, "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCILV("gradSqr3D4n", 3, "shape(-2,-2,-2)", "shape(2,2,2)")

BZ_ET_STENCILM("Jacobian3D", 3, "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCILM("Jacobian3Dn", 3, "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCILM("Jacobian3D4", 3, "shape(-2,-2,-2)", "shape(2,2,2)")
BZ_ET_STENCILM("Jacobian3D4n", 3, "shape(-2,-2,-2)", "shape(2,2,2)")

BZ_ET_STENCIL("curl3D", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCIL("curl3Dn", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCIL("curl3D4", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-2,-2,-2)", "shape(2,2,2)")
BZ_ET_STENCIL("curl3D4n", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-2,-2,-2)", "shape(2,2,2)")
BZ_ET_STENCIL("curl2D", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCIL("curl2Dn", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCIL("curl2D4", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-2,-2)", "shape(2,2)")
BZ_ET_STENCIL("curl2D4n", "T", "typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype", "shape(-2,-2)", "shape(2,2)")

BZ_ET_STENCIL_SCA("div2D", "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCIL_SCA("div2Dn", "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCIL_SCA("div2D4", "shape(-2,-2)", "shape(2,2)")
BZ_ET_STENCIL_SCA("div2D4n", "shape(-2,-2)", "shape(2,2)")
BZ_ET_STENCIL_SCA("div3D", "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCIL_SCA("div3Dn", "shape(-1,-1,-1)", "shape(1,1,1)")
BZ_ET_STENCIL_SCA("div3D4", "shape(-2,-2,-2)", "shape(2,2,2)")


BZ_ET_STENCIL2("div", "double" , "BZ_PROMOTE(typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_numtype)", "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCIL2("divn", "double" , "BZ_PROMOTE(typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_numtype)", "shape(-1,-1)", "shape(1,1)")
BZ_ET_STENCIL2("div4", "double" , "BZ_PROMOTE(typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_numtype)", "shape(-2,-2)", "shape(2,2)")
BZ_ET_STENCIL2("div4n", "double" , "BZ_PROMOTE(typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_numtype)", "shape(-2,-2)", "shape(2,2)")

BZ_ET_STENCIL_DIFF2("mixed22", -1, 1, -1, 1)
BZ_ET_STENCIL_DIFF2("mixed22n", -1, 1, -1, 1)
BZ_ET_STENCIL_DIFF2("mixed24", -2, 2, -2, 2)
BZ_ET_STENCIL_DIFF2("mixed24n", -2, 2, -2, 2)

os.write("""
BZ_NAMESPACE_END
""")
os.close()
