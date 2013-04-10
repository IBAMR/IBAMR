// -*- C++ -*-
/***************************************************************************
 * blitz/array/stencil-et.h  Expression-template-capable stencils
 *
 * $Id: stencil-et.h,v 1.15 2011/03/25 22:41:17 julianc Exp $
 *
 * Copyright (C) 1997-2011 Todd Veldhuizen <tveldhui@acm.org>
 *
 * This file is a part of Blitz.
 *
 * Blitz is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 *
 * Blitz is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Blitz.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Suggestions:          blitz-devel@lists.sourceforge.net
 * Bugs:                 blitz-support@lists.sourceforge.net
 *
 * For more information, please see the Blitz++ Home Page:
 *    https://sourceforge.net/projects/blitz/
 *
 ****************************************************************************/
#ifndef BZ_ARRAY_STENCIL_ET_MACROS_H
#define BZ_ARRAY_STENCIL_ET_MACROS_H

BZ_NAMESPACE(blitz)

/* This file contains the macros that used to declare all
   stencils. They are now declared in the generated file
   stencil-classes.cc. You only need to include this file if you need
   to declare your own stencil ET classes. Note that this file is NOT
   necessary for declaring stencil operators that are applied with
   applyStencil, as described in the documentation. To use these
   macros to declare an ET stencil called "name", you first need to
   declare a stencil operator called "name_stencilop" using
   e.g. BZ_DECLARE_STENCIL_OPERATOR1, and then call the appropriate
   macro here.
*/


/* To avoid matching to the stencil operator in stencilops.h, we must
   explicitly define stencils that operate on arrays. This macro makes
   this slightly less painful for the majority of the stencil classes. */
#define BZ_ET_STENCIL_REDIRECT(name)					\
  template<typename T, int N>						\
   inline _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > \
   name(const Array<T,N>& d1)						\
   { return name(d1.wrap()); }						\
   template<typename T, int N>						\
   inline _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > \
   name(Array<T,N>& d1)							\
   { return name(d1.wrap()); }


/* Defines a stencil ET that operates on an array<P_numtype, N_rank>
   and specifies the return type as array<result, N_rank>. The result
   type is used when running on an array and the etresult type when
   running on an expression. If you want to refer to the native type
   of the expression, set result="P_numtype" and etresult="typename
   T1::T_numtype". Sorry for that ugliness, but they define types
   differently. */

#define BZ_ET_STENCIL(name,result, etresult, MINB, MAXB)		\
  template<typename P_expr, _bz_typename P_numtype>			\
  class name ## _et : public _bz_StencilExpr<P_expr, P_numtype>		\
  {									\
  public:								\
    typedef _bz_StencilExpr<P_expr, P_numtype> T_base;			\
    typedef _bz_typename T_base::T_numtype T_numtype;			\
    typedef _bz_typename T_base::T_expr T_expr;				\
									\
    /* if P_numtype is an ET-type, we need to return an expr */		\
    typedef typename selectET<P_numtype,				\
			      T_numtype,				\
      ETBase<_bz_ArrayExpr<_bz_ArrayExprConstant<P_numtype> > > >::T_selected T_typeprop;\
    typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;	\
    typedef T_numtype T_optype;						\
									\
    /* dummy */								\
    template<int N> struct tvresult {					\
      typedef name ## _et<typename T_expr::template tvresult<N>::Type,T_numtype> Type; \
    };									\
									\
    typedef  name ## _et<_bz_typename P_expr::T_range_result, T_numtype> T_range_result; \
									\
    using T_base::iter_;						\
    using T_base::rank_;						\
  public:								\
    name ## _et(const name ## _et& a) :					\
    _bz_StencilExpr<P_expr, T_numtype>(a)				\
      { }								\
									\
    name ## _et(BZ_ETPARM(T_expr) a) :					\
    _bz_StencilExpr<P_expr, T_numtype>(a)				\
      { }								\
									\
    name ## _et(_bz_typename T_expr::T_ctorArg1 a) :			\
    _bz_StencilExpr<P_expr, T_numtype>(a)				\
      { }								\
									\
    T_result operator*() const						\
    { return name ## _stencilop(iter_); }				\
									\
    /* this is not really const, because we don't undo the moveTo, but	\
       that should not be visible to outside. It would be if we used	\
       some kind of mixed index and stack traversal, but that will	\
       screw things up, const or not. */				\
    template<int N_rank2>						\
      T_result operator()(const TinyVector<int, N_rank2>& i) const	\
    { iter_.moveTo(i); return name ## _stencilop(iter_); }				\
    									\
    T_range_result operator()(const RectDomain<rank_>& d) const		\
    { return T_range_result(iter_(d)); }				\
									\
    T_result operator[](int i) const					\
    { return name ## _stencilop(iter_[i]); }						\
									\
    T_result fastRead(sizeType i) const				\
    {/* this probably isn't very fast... */				\
      iter_._bz_offsetData(i);						\
      T_result r = name ## _stencilop (iter_);					\
      iter_._bz_offsetData(-i);						\
      return r;								\
    }									\
    									\
    /** This way of vectorizing won't work on stencils. */		\
    template<int N>							\
    typename tvresult<N>::Type fastRead_tv(int i) const {		\
      BZPRECHECK(0, "Can't vectorize stencils");			\
      return iter_.template fastRead_tv<N>(i); }					\
									\
    T_result shift(int offset, int dim) const				\
    {									\
      iter_._bz_offsetData(offset, dim);				\
      T_result r = name ## _stencilop (iter_);					\
      iter_._bz_offsetData(-offset, dim);				\
      return r;								\
    }									\
									\
    T_result shift(int offset1, int dim1, int offset2, int dim2) const	\
    {									\
      iter_._bz_offsetData(offset1, dim1, offset2, dim2);		\
      T_result r = name ## _stencilop (iter_);					\
      iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		\
      return r;								\
    }									\
									\
    void prettyPrint(BZ_STD_SCOPE(string) &str,				\
		     prettyPrintFormat& format) const			\
    {									\
      str += "name (stencil)";						\
      str += "(";							\
      iter_.prettyPrint(str, format);					\
      str += ")";							\
    }									\
    									\
    template<typename T1, typename T2 = nilArraySection,		\
      class T3 = nilArraySection, typename T4 = nilArraySection,	\
      class T5 = nilArraySection, typename T6 = nilArraySection,	\
      class T7 = nilArraySection, typename T8 = nilArraySection,	\
      class T9 = nilArraySection, typename T10 = nilArraySection,	\
      class T11 = nilArraySection>					\
      class SliceInfo {							\
      public:								\
      typedef name ## _et<T_expr, T_numtype> T_slice;			\
      };								\
    									\
    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, \
      typename T7, typename T8, typename T9, typename T10, typename T11> \
      name ## _et							\
      operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const \
    {									\
    /* because stencils work inherently in several dimensions it's	\
       complicated to slice the domain. slices will be changed to unit	\
       \ ranges instead. slicing stencil result thus *never* changes	\
       the rank of the expression, unlike the normal case. */		\
      return name ## _et						\
	(iter_(_bz_makeRange(r1),					\
	       _bz_makeRange(r2),					\
	       _bz_makeRange(r3),					\
	       _bz_makeRange(r4),					\
	       _bz_makeRange(r5),					\
	       _bz_makeRange(r6),					\
	       _bz_makeRange(r7),					\
	       _bz_makeRange(r8),					\
	       _bz_makeRange(r9),					\
	       _bz_makeRange(r10),					\
	       _bz_makeRange(r11)));					\
    }									\
									\
  };									\
  /* generate an ET object from an expression */			\
  template<typename T1>							\
  inline _bz_ArrayExpr<name ## _et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, etresult> > \
  name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1)				\
  {									\
    return _bz_ArrayExpr<name ## _et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, etresult> > \
      (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),MINB, MAXB))); \
  }									\
  /* redirect calls with bare arrays to the main function */		\
  template<typename T, int N>						\
  inline _bz_ArrayExpr<name ## _et<typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, result> > \
  name(const Array<T,N>& d1)						\
  { return name(d1.wrap()); }						\
									\
  template<typename T, int N>						\
  inline _bz_ArrayExpr<name ## _et<typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, result> > \
  name(Array<T,N>& d1)							\
  { return name(d1.wrap()); }


/* Defines a stencil ET that operates on two arrays of arbitrary type
   and specifies the return type as array<result, N_rank>. The result
   type is used when running on an array and the etresult type when
   running on an expression. If you want to refer to the native type
   of the expression, set result="P_numtype" and etresult="typename
   T1::T_numtype". Sorry for that ugliness, but they define types
   differently. */

#define BZ_ET_STENCIL2(name,result, etresult, MINB, MAXB)		\
  template<typename P_expr1, typename P_expr2, _bz_typename P_numtype>	\
  class name ## _et2 : public _bz_StencilExpr2<P_expr1, P_expr2, P_numtype> \
  {									\
  public:								\
    typedef _bz_StencilExpr2<P_expr1, P_expr2, P_numtype> T_base;	\
    typedef _bz_typename T_base::T_numtype T_numtype;			\
    typedef _bz_typename T_base::T_expr1 T_expr1;			\
    typedef _bz_typename T_base::T_expr2 T_expr2;			\
									\
    /* if P_numtype is an ET-type, we need to return an expr */		\
    typedef typename selectET<P_numtype,				\
			      T_numtype,				\
      ETBase<_bz_ArrayExpr<_bz_ArrayExprConstant<P_numtype> > > >::T_selected T_typeprop;\
    typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;	\
    typedef T_numtype T_optype;						\
									\
    /* dummy */								\
    template<int N> struct tvresult {					\
      typedef name ## _et2<typename T_expr1::template tvresult<N>::Type,typename T_expr2::template tvresult<N>::Type,T_numtype> Type; \
    };									\
									\
    typedef  name ## _et2<_bz_typename P_expr1::T_range_result, _bz_typename P_expr2::T_range_result, T_numtype> T_range_result; \
									\
    using T_base::iter1_;						\
    using T_base::iter2_;						\
    using T_base::rank_;							\
  public:								\
    name ## _et2(const name ## _et2& a) :				\
      _bz_StencilExpr2<P_expr1, P_expr2, T_numtype>(a)			\
      { }								\
									\
    name ## _et2(BZ_ETPARM(T_expr1) a, BZ_ETPARM(T_expr2) b) :		\
      _bz_StencilExpr2<P_expr1, P_expr2, T_numtype>(a, b)		\
      { }								\
    /*									\
    name ## _et2(_bz_typename T_expr::T_ctorArg1 a) :			\
      _bz_StencilExpr2<P_expr, T_numtype>(a)				\
      { }								\
    */									\
    T_result operator*() const						\
    { return name ## _stencilop(iter1_, iter2_); }					\
									\
    T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const \
    { iter1_.moveTo(i); iter2_.moveTo(i); return name ## _stencilop(iter1_, iter2_); } \
									\
    T_range_result operator()(const RectDomain<rank_>& d) const		\
    { return T_range_result(iter1_(d), iter2_(d)); }			\
									\
    T_result operator[](int i) const					\
    { return name ## _stencilop(iter1_[i], iter2_[i]); }				\
									\
    T_result fastRead(sizeType i) const					\
    {/* this probably isn't very fast... */				\
      iter1_._bz_offsetData(i); iter2_._bz_offsetData(i);		\
      T_result r = name ## _stencilop (iter1_, iter2_);				\
      iter1_._bz_offsetData(-i); iter2_._bz_offsetData(-i);		\
      return r;								\
    }									\
    									\
    /** This way of vectorizing won't work on stencils. */		\
    template<int N>							\
      typename tvresult<N>::Type fastRead_tv(int i) const {		\
      BZPRECHECK(0, "Can't vectorize stencils");			\
      return typename tvresult<N>::Type(iter1_.template fastRead_tv<N>(i),	\
					iter2_.template fastRead_tv<N>(i)); }	\
    									\
    T_result shift(int offset, int dim) const				\
    {									\
      iter1_._bz_offsetData(offset, dim);				\
      iter2_._bz_offsetData(offset, dim);				\
      T_result r = name ## _stencilop (iter1_, iter2_);				\
      iter1_._bz_offsetData(-offset, dim);				\
      iter2_._bz_offsetData(-offset, dim);				\
      return r;								\
    }									\
									\
    T_result shift(int offset1, int dim1, int offset2, int dim2) const	\
    {									\
      iter1_._bz_offsetData(offset1, dim1, offset2, dim2);		\
      iter2_._bz_offsetData(offset1, dim1, offset2, dim2);		\
      T_result r = name ## _stencilop (iter1_, iter2_);				\
      iter1_._bz_offsetData(-offset1, dim1, -offset2, dim2);		\
      iter2_._bz_offsetData(-offset1, dim1, -offset2, dim2);		\
      return r;								\
    }									\
									\
    void prettyPrint(BZ_STD_SCOPE(string) &str,				\
		     prettyPrintFormat& format) const			\
    {									\
      str += "name (stencil)";						\
      str += "(";							\
      iter1_.prettyPrint(str, format);					\
      str += ", ";							\
      iter2_.prettyPrint(str, format);					\
      str += ")";							\
    }									\
									\
    template<typename T1, typename T2 = nilArraySection,		\
      class T3 = nilArraySection, typename T4 = nilArraySection,	\
      class T5 = nilArraySection, typename T6 = nilArraySection,	\
      class T7 = nilArraySection, typename T8 = nilArraySection,	\
      class T9 = nilArraySection, typename T10 = nilArraySection,	\
      class T11 = nilArraySection>					\
      class SliceInfo {							\
      public:								\
      typedef name ## _et2<T_expr1, T_expr2, T_numtype> T_slice;	\
      };								\
    									\
    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, \
      typename T7, typename T8, typename T9, typename T10, typename T11> \
      name ## _et2							\
      operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const \
    {									\
      return name ## _et2						\
	(iter_(_bz_makeRange(r1),					\
	       _bz_makeRange(r2),					\
	       _bz_makeRange(r3),					\
	       _bz_makeRange(r4),					\
	       _bz_makeRange(r5),					\
	       _bz_makeRange(r6),					\
	       _bz_makeRange(r7),					\
	       _bz_makeRange(r8),					\
	       _bz_makeRange(r9),					\
	       _bz_makeRange(r10),					\
	       _bz_makeRange(r11)));					\
    }									\
  };									\
									\
  /* create ET object from application to expression */			\
  template<typename T1, typename T2>					\
  inline _bz_ArrayExpr<name ## _et2<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_range_result, etresult> > \
  name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1,				\
       const BZ_BLITZ_SCOPE(ETBase)<T2>& d2)				\
  {									\
    return _bz_ArrayExpr<name ## _et2<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_range_result, etresult> > \
      (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),MINB, MAXB)), \
       BZ_BLITZ_SCOPE(asExpr)<T2>::getExpr(d2.unwrap())(_bz_shrinkDomain(d2.unwrap().domain(),MINB, MAXB))); \
  }									\
  /* matches to calls involving bare arrays (this is very annoying	\
     because we have to exactly match every possible call combination	\
     to ensure that this matches instead of the operator in		\
     stencilops.h) */							\
  template<typename T1, typename T2, int N2>				\
  inline _bz_ArrayExpr<name ## _et2<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T2,N2> >::T_expr::T_range_result, result> > \
  name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1, Array<T2,N2>& d2)		\
  { return name(d1.wrap(), d2.wrap()); }				\
									\
  template<typename T1, typename T2, int N2>				\
  inline _bz_ArrayExpr<name ## _et2<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T2,N2> >::T_expr::T_range_result, result> > \
  name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1, const Array<T2,N2>& d2)	\
  { return name(d1.wrap(), d2.wrap()); }				\
									\
  template<typename T1, int N1, typename T2>				\
  inline _bz_ArrayExpr<name ## _et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T1,N1> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_range_result, result> > \
  name(Array<T1,N1>& d1, const BZ_BLITZ_SCOPE(ETBase)<T2>& d2)		\
  { return name(d1.wrap(), d2.wrap()); }				\
									\
  template<typename T1, int N1, typename T2>				\
  inline _bz_ArrayExpr<name ## _et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T1,N1> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_range_result, result> > \
  name(const Array<T1,N1>& d1, const BZ_BLITZ_SCOPE(ETBase)<T2>& d2)	\
  { return name(d1.wrap(), d2.wrap()); }				\
									\
  template<typename T1, int N1, typename T2, int N2>			\
  inline _bz_ArrayExpr<name ## _et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T1,N1> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T2,N2> >::T_expr::T_range_result, result> > \
  name(const Array<T1,N1>& d1, Array<T2,N2>& d2)				\
  { return name(d1.wrap(), d2.wrap()); }				\
									\
  template<typename T1, int N1, typename T2, int N2>			\
  inline _bz_ArrayExpr<name ## _et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T1,N1> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T2,N2> >::T_expr::T_range_result, result> > \
  name(Array<T1,N1>& d1, const Array<T2,N2>& d2)			\
  { return name(d1.wrap(), d2.wrap()); }				\
  									\
  template<typename T, int N>						\
  inline _bz_ArrayExpr<name ## _et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, result> > \
  name(Array<T,N>& d1, Array<T,N>& d2)					\
  { return name(d1.wrap(), d2.wrap()); }				\
									\
  template<typename T, int N>						\
  inline _bz_ArrayExpr<name ## _et2<typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result, result> > \
  name(const Array<T,N>& d1, const Array<T,N>& d2)			\
  { return name(d1.wrap(), d2.wrap()); }


 /* Defines a stencil ET that operates on an array<P_numtype, N_rank>
    and returns a multicomponent array<TinyMatrix<P_numtype::T_element,
    rank, rank> >, N_rank>. P_numtype can be a TinyVector or a scalar,
    I think. */

 #define BZ_ET_STENCILM(name,result_rank, MINB, MAXB)			\
   template<typename P_expr>						\
   class name ## _et : public _bz_StencilExpr<P_expr, TinyMatrix<_bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element, result_rank, result_rank> > \
   {									\
 public:									\
     typedef _bz_StencilExpr<P_expr, TinyMatrix<_bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element, result_rank, result_rank> > T_base; \
     typedef _bz_typename T_base::T_numtype T_numtype;			\
     typedef _bz_typename T_base::T_expr T_expr;			\
     									\
     /* there is no return type selection, as we are returning a	\
	TinyMatrix. This must be returned as a FastTMCopyIterator since the \
	output of the stencil operator is a temporary. */		\
  typedef ETBase<_bz_ArrayExpr<FastTM2CopyIterator<typename multicomponent_traits<typename P_expr::T_numtype>::T_element, result_rank, result_rank> > > T_typeprop;\
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;\
  typedef T_numtype T_optype;\
									\
    /* dummy */								\
    template<int N> struct tvresult {					\
      typedef name ## _et<typename T_expr::template tvresult<N>::Type> Type; \
    };									\
  									\
     typedef  name ## _et<_bz_typename P_expr::T_range_result> T_range_result; \
									 \
     using T_base::iter_;						\
     using T_base::rank_;						\
   public:								\
     name ## _et(const name ## _et& a) :				\
       _bz_StencilExpr<P_expr, T_numtype>(a)				\
       { }								\
     									\
     name ## _et(BZ_ETPARM(T_expr) a) :					\
       _bz_StencilExpr<P_expr, T_numtype>(a)				\
       { }								\
									 \
     name ## _et(_bz_typename T_expr::T_ctorArg1 a) :			\
       _bz_StencilExpr<P_expr, T_numtype>(a)				\
       { }								\
     									\
     T_result operator*() const					\
     { return name ## _stencilop(iter_); }						\
     T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const \
     { iter_.moveTo(i); return name ## _stencilop(iter_); }				\
     									\
     T_range_result operator()(const RectDomain<rank_>& d) const		\
     { return T_range_result(iter_(d)); }				\
     									\
     T_result operator[](int i) const					\
     { return name ## _stencilop(iter_[i]); }						\
     									\
     T_result fastRead(sizeType i) const				\
     {/* this probably isn't very fast... */				\
       iter_._bz_offsetData(i);						\
       T_result r = name ## _stencilop (iter_);					\
       iter_._bz_offsetData(-i);					\
       return r;							\
     }									\
    									\
    /** This way of vectorizing won't work on stencils. */		\
    template<int N>							\
    typename tvresult<N>::Type fastRead_tv(int i) const {		\
      BZPRECHECK(0, "Can't vectorize stencils");			\
      return iter_.template fastRead_tv<N>(i); }					\
									\
     T_result shift(int offset, int dim) const				\
     {									\
       iter_._bz_offsetData(offset, dim);				\
       T_result r = name ## _stencilop (iter_);					\
       iter_._bz_offsetData(-offset, dim);				\
       return r;							\
     }									\
									 \
     T_result shift(int offset1, int dim1, int offset2, int dim2) const \
     {									\
       iter_._bz_offsetData(offset1, dim1, offset2, dim2);		\
       T_result r = name ## _stencilop (iter_);					\
       iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		\
       return r;							\
     }									\
									 \
     void prettyPrint(BZ_STD_SCOPE(string) &str,			\
		      prettyPrintFormat& format) const			\
     {									\
       str += "name (stencil)";						\
       str += "(";							\
       iter_.prettyPrint(str, format);					\
       str += ")";							\
     }									\
									\
    template<typename T1, typename T2 = nilArraySection,		\
      class T3 = nilArraySection, typename T4 = nilArraySection,	\
      class T5 = nilArraySection, typename T6 = nilArraySection,	\
      class T7 = nilArraySection, typename T8 = nilArraySection,	\
      class T9 = nilArraySection, typename T10 = nilArraySection,	\
      class T11 = nilArraySection>					\
      class SliceInfo {							\
      public:								\
      typedef name ## _et<T_expr> T_slice;				\
      };								\
    									\
    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, \
      typename T7, typename T8, typename T9, typename T10, typename T11> \
      name ## _et							\
      operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const \
    {									\
      return name ## _et						\
	(iter_(_bz_makeRange(r1),					\
	       _bz_makeRange(r2),					\
	       _bz_makeRange(r3),					\
	       _bz_makeRange(r4),					\
	       _bz_makeRange(r5),					\
	       _bz_makeRange(r6),					\
	       _bz_makeRange(r7),					\
	       _bz_makeRange(r8),					\
	       _bz_makeRange(r9),					\
	       _bz_makeRange(r10),					\
	       _bz_makeRange(r11)));					\
    }									\
   };									\
   /* create ET from application to expression */			\
   template<typename T1>						\
   inline _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> > \
   name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1)				\
   {									\
     return _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	\
       (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),MINB, MAXB))); \
       }								\
   BZ_ET_STENCIL_REDIRECT(name)						\


 /* Defines a stencil ET that operates on a (scalar) array<P_numtype,
    N_rank> and returns a multicomponent
    array<TinyVector<P_numtype::T_element, result_rank> >, N_rank>. */

 #define BZ_ET_STENCILV(name,result_rank, MINB, MAXB)			\
   template<typename P_expr>						\
   class name ## _et : public _bz_StencilExpr<P_expr, TinyVector<typename P_expr::T_numtype,result_rank> > \
   {									\
 public:									\
     typedef _bz_StencilExpr<P_expr, TinyVector<typename P_expr::T_numtype,result_rank> > T_base; \
     typedef _bz_typename T_base::T_numtype T_numtype;			\
     typedef _bz_typename T_base::T_expr T_expr;			\
									\
     /* there is no return type selection, we assume P_numtype is scalar \
	and that we are returning a TinyVector. This needs to be returned \
	as a FastTVCopyIterator that keeps a copy of the TV it is	\
	iterating over, since the result of the stencil operator is	\
	a temporary. */							\
     typedef ETBase<_bz_ArrayExpr<FastTV2CopyIterator<typename P_expr::T_numtype, result_rank> > > T_typeprop; \
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;		\
  typedef typename T_expr::T_numtype T_optype;				\
									\
    /* dummy */								\
    template<int N> struct tvresult {					\
      typedef name ## _et<typename T_expr::template tvresult<N>::Type> Type; \
    };									\
									\
     typedef  name ## _et<_bz_typename P_expr::T_range_result> T_range_result; \
									 \
     using T_base::iter_;						\
     using T_base::rank_;							\
   public:								\
     name ## _et(const name ## _et& a) :					\
     _bz_StencilExpr<P_expr, T_numtype>(a)				\
       { }								\
									 \
     name ## _et(BZ_ETPARM(T_expr) a) :					\
     _bz_StencilExpr<P_expr, T_numtype>(a)				\
       { }								\
									 \
     name ## _et(_bz_typename T_expr::T_ctorArg1 a) :			\
     _bz_StencilExpr<P_expr, T_numtype>(a)				\
       { }								\
									 \
     T_result operator*() const					\
     { return name ## _stencilop(iter_); }						\
     T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const \
     { iter_.moveTo(i); return name ## _stencilop(iter_); }				\
									 \
     T_range_result operator()(const RectDomain<rank_>& d) const		\
     { return T_range_result(iter_(d)); }				\
									 \
     T_result operator[](int i) const					\
     { return name ## _stencilop(iter_[i]); }						\
									 \
     T_result fastRead(sizeType i) const				\
     {/* this probably isn't very fast... */				\
       iter_._bz_offsetData(i);						\
       T_result r = name ## _stencilop (iter_);					\
       iter_._bz_offsetData(-i);					\
       return r;							\
     }									\
    									\
    /** This way of vectorizing won't work on stencils. */		\
    template<int N>							\
    typename tvresult<N>::Type fastRead_tv(int i) const {		\
      BZPRECHECK(0, "Can't vectorize stencils");			\
      return iter_.template fastRead_tv<N>(i); }					\
									 \
     T_result shift(int offset, int dim) const				\
     {									\
       iter_._bz_offsetData(offset, dim);				\
       T_result r = name ## _stencilop (iter_);					\
       iter_._bz_offsetData(-offset, dim);				\
       return r;								\
     }									\
									 \
     T_result shift(int offset1, int dim1, int offset2, int dim2) const \
     {									\
       iter_._bz_offsetData(offset1, dim1, offset2, dim2);		\
       T_result r = name ## _stencilop (iter_);					\
       iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		\
       return r;								\
     }									\
									 \
     void prettyPrint(BZ_STD_SCOPE(string) &str,				\
		      prettyPrintFormat& format) const			\
     {									\
       str += "name (stencil)";						\
       str += "(";							\
       iter_.prettyPrint(str, format);					\
       str += ")";							\
     }									\
									\
    template<typename T1, typename T2 = nilArraySection,		\
      class T3 = nilArraySection, typename T4 = nilArraySection,	\
      class T5 = nilArraySection, typename T6 = nilArraySection,	\
      class T7 = nilArraySection, typename T8 = nilArraySection,	\
      class T9 = nilArraySection, typename T10 = nilArraySection,	\
      class T11 = nilArraySection>					\
      class SliceInfo {							\
      public:								\
      typedef name ## _et<T_expr> T_slice;				\
      };								\
    									\
    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, \
      typename T7, typename T8, typename T9, typename T10, typename T11> \
      name ## _et							\
      operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const \
    {									\
      return name ## _et						\
	(iter_(_bz_makeRange(r1),					\
	       _bz_makeRange(r2),					\
	       _bz_makeRange(r3),					\
	       _bz_makeRange(r4),					\
	       _bz_makeRange(r5),					\
	       _bz_makeRange(r6),					\
	       _bz_makeRange(r7),					\
	       _bz_makeRange(r8),					\
	       _bz_makeRange(r9),					\
	       _bz_makeRange(r10),					\
	       _bz_makeRange(r11)));					\
    }									\
   };									\
  /* create ET from application to expression */				\
   template<typename T1>							\
   inline _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	\
   name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1)				\
   {									\
     return _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	\
       (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),MINB, MAXB))); \
   }									\
   BZ_ET_STENCIL_REDIRECT(name)


 /* Defines a stencil ET that operates on an array<P_numtype, N_rank>
    (where P_numtype presumably is a multicomponent type) and returns a
    scalar array<P_numtype::T_element, N_rank>. */

 #define BZ_ET_STENCIL_SCA(name, MINB, MAXB)				\
   template<typename P_expr>						\
   class name ## _et : public _bz_StencilExpr<P_expr, _bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element> \
   {									\
 public:									\
     typedef _bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element T_result; \
     typedef _bz_StencilExpr<P_expr, T_result> T_base;			\
     typedef _bz_typename T_base::T_numtype T_numtype;			\
     typedef _bz_typename T_base::T_expr T_expr;			\
     									\
     /* there is no selecting return type here. because we *know* it is	\
	scalar T_result, there's no question of whether we could be doing \
	multicomponent evaluations. */					\
     typedef T_result T_typeprop;					\
     typedef T_numtype T_optype;					\
									\
    /* dummy */								\
    template<int N> struct tvresult {					\
      typedef name ## _et<typename T_expr::template tvresult<N>::Type> Type; \
    };									\
									\
     typedef  name ## _et<_bz_typename P_expr::T_range_result> T_range_result; \
     									\
     using T_base::iter_;						\
     using T_base::rank_;							\
   public:								\
     name ## _et(const name ## _et& a) :					\
     _bz_StencilExpr<P_expr, T_numtype>(a)				\
       { }								\
									 \
     name ## _et(BZ_ETPARM(T_expr) a) :					\
     _bz_StencilExpr<P_expr, T_numtype>(a)				\
       { }								\
									 \
     name ## _et(_bz_typename T_expr::T_ctorArg1 a) :			\
     _bz_StencilExpr<P_expr, T_numtype>(a)				\
             { }								\
    									\
    T_result operator*() const						\
    { return name ## _stencilop(iter_); }						\
    T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const \
    { iter_.moveTo(i); return name ## _stencilop(iter_); }				\
									\
    T_range_result operator()(const RectDomain<rank_>& d) const		\
    { return T_range_result(iter_(d)); }				\
									\
    T_result operator[](int i) const					\
    { return name ## _stencilop(iter_[i]); }						\
									\
    T_result fastRead(sizeType i) const				\
    {/* this probably isn't very fast... */				\
      iter_._bz_offsetData(i);						\
      T_result r = name ## _stencilop (iter_);					\
      iter_._bz_offsetData(-i);						\
      return r;								\
    }									\
    									\
    /** This way of vectorizing won't work on stencils. */		\
    template<int N>							\
    typename tvresult<N>::Type fastRead_tv(int i) const {		\
      BZPRECHECK(0, "Can't vectorize stencils");			\
      return iter_.template fastRead_tv<N>(i); }					\
									\
    T_result shift(int offset, int dim) const				\
    {									\
      iter_._bz_offsetData(offset, dim);				\
      T_result r = name ## _stencilop (iter_);					\
      iter_._bz_offsetData(-offset, dim);				\
      return r;								\
    }									\
									\
    T_result shift(int offset1, int dim1, int offset2, int dim2) const	\
    {									\
      iter_._bz_offsetData(offset1, dim1, offset2, dim2);		\
      T_result r = name ## _stencilop (iter_);					\
      iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		\
      return r;								\
    }									\
									\
    void prettyPrint(BZ_STD_SCOPE(string) &str,				\
		     prettyPrintFormat& format) const			\
    {									\
      str += "name (stencil)";						\
      str += "(";							\
      iter_.prettyPrint(str, format);					\
      str += ")";							\
    }									\
									\
    template<typename T1, typename T2 = nilArraySection,		\
      class T3 = nilArraySection, typename T4 = nilArraySection,	\
      class T5 = nilArraySection, typename T6 = nilArraySection,	\
      class T7 = nilArraySection, typename T8 = nilArraySection,	\
      class T9 = nilArraySection, typename T10 = nilArraySection,	\
      class T11 = nilArraySection>					\
      class SliceInfo {							\
      public:								\
      typedef name ## _et<T_expr> T_slice;				\
      };								\
    									\
    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, \
      typename T7, typename T8, typename T9, typename T10, typename T11> \
      name ## _et							\
      operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const \
    {									\
      return name ## _et						\
	(iter_(_bz_makeRange(r1),					\
	       _bz_makeRange(r2),					\
	       _bz_makeRange(r3),					\
	       _bz_makeRange(r4),					\
	       _bz_makeRange(r5),					\
	       _bz_makeRange(r6),					\
	       _bz_makeRange(r7),					\
	       _bz_makeRange(r8),					\
	       _bz_makeRange(r9),					\
	       _bz_makeRange(r10),					\
	       _bz_makeRange(r11)));					\
    }									\
  };									\
  /* create ET from application to expression */			\
  template<typename T1>							\
  inline _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	\
  name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1)				\
  {									\
    return _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	\
      (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),MINB, MAXB))); \
  }									\
   BZ_ET_STENCIL_REDIRECT(name)


/* Defines a stencil ET difference operator that operates on an
   array<P_numtype, N_rank> and returns a array<P_numtype,
   N_rank>. (The only significance of the "difference" aspect is that
   the operator is assumed to take a second argument which is the
   dimension to do the difference in). MINB and MAXB are integer
   expressions describing the extent of the operator in the operating
   dimension. */

#define BZ_ET_STENCIL_DIFF(name, MINB, MAXB)				\
  template<typename P_expr>						\
  class name ## _et :							\
    public _bz_StencilExpr<P_expr, _bz_typename P_expr::T_numtype>	\
  {									\
  public:								\
    typedef _bz_StencilExpr<P_expr, _bz_typename P_expr::T_numtype> T_base; \
    typedef _bz_typename T_base::T_numtype T_numtype;			\
    typedef _bz_typename T_base::T_expr T_expr;				\
    									\
    /* select return type */						\
    typedef typename unwrapET<typename T_expr::T_result>::T_unwrapped test; \
    typedef typename selectET<typename T_expr::T_typeprop,		\
			      T_numtype,				\
			      name ##_et<test> >::T_selected T_typeprop; \
    typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;	\
    typedef T_numtype T_optype;						\
									\
    /* dummy */								\
    template<int N> struct tvresult {					\
      typedef name ## _et<typename T_expr::template tvresult<N>::Type> Type; \
    };									\
    									\
    typedef  name ## _et<_bz_typename P_expr::T_range_result> T_range_result; \
    									\
    using T_base::iter_;						\
    using T_base::rank_;						\
  public:								\
    name ## _et(const name ## _et& a) :					\
      _bz_StencilExpr<P_expr, T_numtype>(a), dim_(a.dim_)		\
      { }								\
    									\
    name ## _et(BZ_ETPARM(T_expr) a, int dim) :				\
      _bz_StencilExpr<P_expr, T_numtype>(a), dim_(dim)			\
    { }									\
    									\
    name ## _et(_bz_typename T_expr::T_ctorArg1 a, int dim) :		\
      _bz_StencilExpr<P_expr, T_numtype>(a), dim_(dim)			\
    { }									\
    									\
    T_result operator*() const						\
    { return name ## _stencilop(iter_, dim_); }				\
    T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const \
    { iter_.moveTo(i); return name ## _stencilop(iter_, dim_); }	\
    									\
    T_range_result operator()(const RectDomain<rank_>& d) const		\
    { return T_range_result(iter_(d), dim_); }				\
    									\
    T_result operator[](int i) const					\
    { return name ## _stencilop(iter_[i], dim_); }			\
    									\
    T_result fastRead(sizeType i) const					\
    {/* this probably isn't very fast... */				\
      iter_._bz_offsetData(i);						\
      T_result r = name ## _stencilop (iter_, dim_);			\
      iter_._bz_offsetData(-i);						\
      return r;								\
    }									\
    									\
    /** This way of vectorizing won't work on stencils. */		\
    template<int N>							\
    typename tvresult<N>::Type fastRead_tv(int i) const {		\
      BZPRECHECK(0, "Can't vectorize stencils");			\
      return typename tvresult<N>::Type(iter_.template fastRead_tv<N>(i),dim_); } \
									\
    T_result shift(int offset, int dim) const				\
    {									\
      iter_._bz_offsetData(offset, dim);				\
      T_result r = name ## _stencilop (iter_);				\
      iter_._bz_offsetData(-offset, dim);				\
      return r;								\
    }									\
    									\
    T_result shift(int offset1, int dim1, int offset2, int dim2) const	\
    {									\
      iter_._bz_offsetData(offset1, dim1, offset2, dim2);		\
      T_result r = name ## _stencilop (iter_);				\
      iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		\
      return r;								\
    }									\
									\
    void prettyPrint(BZ_STD_SCOPE(string) &str,				\
		     prettyPrintFormat& format) const			\
    {									\
      str += "name (stencil)";						\
      str += "(";							\
      iter_.prettyPrint(str, format);					\
      str += ")";							\
    }									\
									\
    template<typename T1, typename T2 = nilArraySection,		\
      class T3 = nilArraySection, typename T4 = nilArraySection,	\
      class T5 = nilArraySection, typename T6 = nilArraySection,	\
      class T7 = nilArraySection, typename T8 = nilArraySection,	\
      class T9 = nilArraySection, typename T10 = nilArraySection,	\
      class T11 = nilArraySection>					\
      class SliceInfo {							\
      public:								\
      typedef name ## _et<T_expr> T_slice;				\
      };								\
    									\
    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, \
      typename T7, typename T8, typename T9, typename T10, typename T11> \
      name ## _et							\
      operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const \
    {									\
      return name ## _et						\
	(iter_(_bz_makeRange(r1),					\
	       _bz_makeRange(r2),					\
	       _bz_makeRange(r3),					\
	       _bz_makeRange(r4),					\
	       _bz_makeRange(r5),					\
	       _bz_makeRange(r6),					\
	       _bz_makeRange(r7),					\
	       _bz_makeRange(r8),					\
	       _bz_makeRange(r9),					\
	       _bz_makeRange(r10),					\
	       _bz_makeRange(r11)),dim_);				\
    }									\
									\
  private:								\
    int dim_;								\
  };									\
 /* create ET from application to expression */				\
  template<typename T1>							\
  inline _bz_ArrayExpr<name ## _et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	\
  name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1, int dim)			\
  {									\
    TinyVector<int, BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::rank_> minb(0), maxb(0); \
    minb[dim]=MINB; maxb[dim]=MAXB;					\
    return _bz_ArrayExpr<name ## _et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	\
      (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),minb, maxb)), dim); \
  }									\
  /* forward operations on arrays to main function */			\
  template<typename T, int N>						\
   inline _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > \
  name(const Array<T,N>& d1, int dim)					\
  { return name(d1.wrap(), dim); }					\
   template<typename T, int N>						\
   inline _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > \
   name(Array<T,N>& d1, int dim)					\
   { return name(d1.wrap(), dim); }


/* Defines a stencil ET difference operator that operates on a
   multicomponent array<P_numtype, N_rank> and returns an
   array<P_numtype::T_element, N_rank>. */

#define BZ_ET_STENCIL_MULTIDIFF(name, MINB, MAXB)			\
  template<typename P_expr>						\
  class name ## _et_multi : public _bz_StencilExpr<P_expr, _bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element> \
  {									\
  public:								\
    typedef _bz_typename multicomponent_traits<typename P_expr::T_numtype>::T_element T_result; \
    typedef _bz_StencilExpr<P_expr, T_result> T_base;		\
    typedef _bz_typename T_base::T_numtype T_numtype;			\
    typedef _bz_typename T_base::T_expr T_expr;				\
									\
    /* there is no selecting return type here. because we *know* it is	\
       T_result, there's no question of whether we could be doing	\
       multicomponent evaluations. */					\
    typedef T_result T_typeprop;					\
  typedef T_numtype T_optype;						\
									\
    /* dummy */								\
    template<int N> struct tvresult {					\
      typedef name ## _et_multi<typename T_expr::template tvresult<N>::Type> Type; \
    };									\
									\
    typedef  name ## _et_multi<_bz_typename P_expr::T_range_result> T_range_result; \
									\
    using T_base::iter_;						\
    using T_base::rank_;							\
  public:								\
    name ## _et_multi(const name ## _et_multi& a) :		\
      _bz_StencilExpr<P_expr, T_numtype>(a), comp_(a.comp_), dim_(a.dim_) \
      { }								\
									\
    name ## _et_multi(BZ_ETPARM(T_expr) a, int comp, int dim) :	\
      _bz_StencilExpr<P_expr, T_numtype>(a),			\
      comp_(comp), dim_(dim)						\
      { }								\
									\
    name ## _et_multi(_bz_typename T_expr::T_ctorArg1 a, int comp, int dim) : \
      _bz_StencilExpr<P_expr, T_numtype>(a),			\
      comp_(comp), dim_(dim)						\
      { }								\
									\
    T_result operator*() const						\
    { return name ## _stencilop(iter_, comp_, dim_); }				\
    T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const \
    { iter_.moveTo(i); return name ## _stencilop(iter_, comp_, dim_); }		\
									\
    T_range_result operator()(const RectDomain<rank_>& d) const		\
    { return T_range_result(iter_(d), comp_, dim_); }			\
									\
    T_result operator[](int i) const					\
    { return name ## _stencilop(iter_[i], comp_, dim_); }				\
									\
    T_result fastRead(sizeType i) const				\
    {/* this probably isn't very fast... */				\
      iter_._bz_offsetData(i);						\
      T_result r = name ## _stencilop (iter_, comp_, dim_);				\
      iter_._bz_offsetData(-i);						\
      return r;								\
    }									\
    									\
    /** This way of vectorizing won't work on stencils. */		\
    template<int N>							\
    typename tvresult<N>::Type fastRead_tv(int i) const {		\
      BZPRECHECK(0, "Can't vectorize stencils");			\
      return typename tvresult<N>::Type(iter_.template fastRead_tv<N>(i),comp_,dim_); } \
									\
    T_result shift(int offset, int dim) const				\
    {									\
      iter_._bz_offsetData(offset, dim);				\
      T_result r = name ## _stencilop (iter_);					\
      iter_._bz_offsetData(-offset, dim);				\
      return r;								\
    }									\
									\
    T_result shift(int offset1, int dim1, int offset2, int dim2) const	\
    {									\
      iter_._bz_offsetData(offset1, dim1, offset2, dim2);		\
      T_result r = name ## _stencilop (iter_);					\
      iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		\
      return r;								\
    }									\
									\
    void prettyPrint(BZ_STD_SCOPE(string) &str,				\
		     prettyPrintFormat& format) const			\
    {									\
      str += "name (stencil)";						\
      str += "(";							\
      iter_.prettyPrint(str, format);					\
      str += ")";							\
    }									\
									\
    template<typename T1, typename T2 = nilArraySection,		\
      class T3 = nilArraySection, typename T4 = nilArraySection,	\
      class T5 = nilArraySection, typename T6 = nilArraySection,	\
      class T7 = nilArraySection, typename T8 = nilArraySection,	\
      class T9 = nilArraySection, typename T10 = nilArraySection,	\
      class T11 = nilArraySection>					\
      class SliceInfo {							\
      public:								\
      typedef name ## _et_multi<T_expr> T_slice;				\
      };								\
    									\
    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, \
      typename T7, typename T8, typename T9, typename T10, typename T11> \
      name ## _et_multi							\
      operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const \
    {									\
      return name ## _et_multi						\
	(iter_(_bz_makeRange(r1),					\
	       _bz_makeRange(r2),					\
	       _bz_makeRange(r3),					\
	       _bz_makeRange(r4),					\
	       _bz_makeRange(r5),					\
	       _bz_makeRange(r6),					\
	       _bz_makeRange(r7),					\
	       _bz_makeRange(r8),					\
	       _bz_makeRange(r9),					\
	       _bz_makeRange(r10),					\
	       _bz_makeRange(r11)),comp_, dim_);			\
    }									\
									\
  private:								\
    int comp_;								\
    int dim_;								\
    };									\
  /* create ET from application to expression */			\
  template<typename T1>							\
  inline _bz_ArrayExpr<name ## _et_multi<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	\
  name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1, int comp, int dim)		\
  {									\
    TinyVector<int, BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::rank_> minb(0), maxb(0); \
    minb[dim]=MINB; maxb[dim]=MAXB;					\
    return _bz_ArrayExpr<name ## _et_multi<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> > \
      (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),minb, maxb)), comp, dim); \
  }									\
  /* forward operations on arrays to main function */			\
  template<typename T, int N>						\
  inline _bz_ArrayExpr<name ## _et_multi<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > \
  name(const Array<T,N>& d1, int comp, int dim)				\
  { return name(d1.wrap(), comp, dim); }				\
									\
  template<typename T, int N>						\
  inline _bz_ArrayExpr<name ## _et_multi<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > \
  name(Array<T,N>& d1, int comp, int dim)				\
  { return name(d1.wrap(), comp, dim); }


/* Defines a stencil ET double-difference operator that operates on an
   array<P_numtype, N_rank> and returns a array<P_numtype,
   N_rank>. (The only significance of the "difference" aspect is that
   the operator is assumed to take two extra arguments which are the
   dimensions to do the differences in). */

#define BZ_ET_STENCIL_DIFF2(name, MINB1, MAXB1, MINB2, MAXB2)		\
 template<typename P_expr>						\
 class name ## _et : public _bz_StencilExpr<P_expr, _bz_typename P_expr::T_numtype> \
 {									\
 public:								\
   typedef _bz_StencilExpr<P_expr, _bz_typename P_expr::T_numtype> T_base;	\
   typedef _bz_typename T_base::T_numtype T_numtype;			\
   typedef _bz_typename T_base::T_expr T_expr;				\
									\
   /* select return type */						\
   typedef typename unwrapET<typename T_expr::T_result>::T_unwrapped test; \
  typedef typename selectET<typename T_expr::T_typeprop,		\
			    T_numtype,					\
			    name ## _et<test> >::T_selected T_typeprop;	\
  typedef typename unwrapET<T_typeprop>::T_unwrapped T_result;		\
  typedef T_numtype T_optype;						\
									\
    /* dummy */								\
    template<int N> struct tvresult {					\
      typedef name ## _et<typename T_expr::template tvresult<N>::Type> Type; \
    };									\
									\
   typedef  name ## _et<_bz_typename P_expr::T_range_result> T_range_result; \
   									\
   using T_base::iter_;							\
   using T_base::rank_;							\
 public:								\
   name ## _et(const name ## _et& a) :					\
   _bz_StencilExpr<P_expr, T_numtype>(a),				\
     dim1_(a.dim1_), dim2_(a.dim2_)					\
     { }								\
   									\
   name ## _et(BZ_ETPARM(T_expr) a, int dim1, int dim2) :		\
   _bz_StencilExpr<P_expr, T_numtype>(a),				\
     dim1_(dim1), dim2_(dim2)						\
   { }									\
   									\
   name ## _et(_bz_typename T_expr::T_ctorArg1 a,			\
	       int dim1, int dim2) :					\
   _bz_StencilExpr<P_expr, T_numtype>(a),				\
     dim1_(dim1), dim2_(dim2)						\
   { }									\
   									\
   T_result operator*() const						\
   { return name ## _stencilop(iter_, dim1_, dim2_); }				\
   T_result operator()(_bz_typename _bz_IndexParameter<TinyVector<int, rank_> >::type i) const \
   { iter_.moveTo(i); return name ## _stencilop(iter_, dim1_, dim2_); }		\
									\
   T_range_result operator()(const RectDomain<rank_>& d) const		\
   { return T_range_result(iter_(d), dim1_, dim2_); }			\
   									\
   T_result operator[](int i) const					\
   { return name ## _stencilop(iter_[i], dim1_, dim2_); }				\
									\
   T_result fastRead(sizeType i) const					\
   {/* this probably isn't very fast... */				\
     iter_._bz_offsetData(i);						\
     T_result r = name ## _stencilop (iter_, dim1_, dim2_);				\
     iter_._bz_offsetData(-i);						\
     return r;								\
   }									\
    									\
    /** This way of vectorizing won't work on stencils. */		\
    template<int N>							\
    typename tvresult<N>::Type fastRead_tv(int i) const {		\
      BZPRECHECK(0, "Can't vectorize stencils");			\
      return typename tvresult<N>::Type(iter_.template fastRead_tv<N>(i),dim1_,dim2_); } \
									\
    T_result shift(int offset, int dim) const				\
    {									\
      iter_._bz_offsetData(offset, dim);				\
      T_result r = name ## _stencilop (iter_);					\
      iter_._bz_offsetData(-offset, dim);				\
      return r;								\
    }									\
									\
    T_result shift(int offset1, int dim1, int offset2, int dim2) const	\
    {									\
      iter_._bz_offsetData(offset1, dim1, offset2, dim2);		\
      T_result r = name ## _stencilop (iter_);					\
      iter_._bz_offsetData(-offset1, dim1, -offset2, dim2);		\
      return r;								\
    }									\
									\
    void prettyPrint(BZ_STD_SCOPE(string) &str,				\
		     prettyPrintFormat& format) const			\
    {									\
      str += "name (stencil)";						\
      str += "(";							\
      iter_.prettyPrint(str, format);					\
      str += ")";							\
    }									\
									\
    template<typename T1, typename T2 = nilArraySection,		\
      class T3 = nilArraySection, typename T4 = nilArraySection,	\
      class T5 = nilArraySection, typename T6 = nilArraySection,	\
      class T7 = nilArraySection, typename T8 = nilArraySection,	\
      class T9 = nilArraySection, typename T10 = nilArraySection,	\
      class T11 = nilArraySection>					\
      class SliceInfo {							\
      public:								\
      typedef name ## _et<T_expr> T_slice;				\
      };								\
    									\
    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, \
      typename T7, typename T8, typename T9, typename T10, typename T11> \
      name ## _et							\
      operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6, T7 r7, T8 r8, T9 r9, T10 r10, T11 r11) const \
    {									\
      return name ## _et						\
	(iter_(_bz_makeRange(r1),					\
	       _bz_makeRange(r2),					\
	       _bz_makeRange(r3),					\
	       _bz_makeRange(r4),					\
	       _bz_makeRange(r5),					\
	       _bz_makeRange(r6),					\
	       _bz_makeRange(r7),					\
	       _bz_makeRange(r8),					\
	       _bz_makeRange(r9),					\
	       _bz_makeRange(r10),					\
	       _bz_makeRange(r11)), dim1_, dim2_);			\
    }									\
   									\
private:								\
   int dim1_, dim2_;							\
 };									\
 									\
 /* create ET from application to expression */				\
template<typename T1>							\
 inline _bz_ArrayExpr<name ## _et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	\
 name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1, int dim1, int dim2)		\
 {									\
   TinyVector<int, BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::rank_> minb(0), maxb(0);	\
   minb[dim1]=MINB1; maxb[dim1]=MAXB1;					\
   minb[dim2]=MINB2; maxb[dim2]=MAXB2;					\
   return _bz_ArrayExpr<name ## _et<typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_range_result> >	\
     (BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap())(_bz_shrinkDomain(d1.unwrap().domain(),minb, maxb)), dim1, dim2); \
 }									\
  /* forward operations on arrays to main function */			\
  template<typename T, int N>						\
  inline _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > \
  name(const Array<T,N>& d1, int dim1, int dim2)			\
  { return name(d1.wrap(), dim1, dim2); }				\
									\
  template<typename T, int N>						\
  inline _bz_ArrayExpr<name ## _et<_bz_typename BZ_BLITZ_SCOPE(asExpr)<Array<T,N> >::T_expr::T_range_result> > \
  name(Array<T,N>& d1, int dim1, int dim2)				\
  { return name(d1.wrap(), dim1, dim2); }

BZ_NAMESPACE_END

#endif // BZ_ARRAY_STENCIL_ET_MACROS_H
