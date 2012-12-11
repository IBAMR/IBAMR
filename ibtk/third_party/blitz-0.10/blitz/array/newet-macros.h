// -*- C++ -*-
/***************************************************************************
 * blitz/array/newet-macros.h  Macros for new e.t. implementation
 *
 * $Id$
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
#ifndef BZ_NEWET_MACROS_H
#define BZ_NEWET_MACROS_H

#include <blitz/array/asexpr.h>

BZ_NAMESPACE(blitz)

#ifdef BZ_HAVE_TEMPLATES_AS_TEMPLATE_ARGUMENTS

/*
 * Unary functions and operators
 */

#define BZ_DECLARE_ARRAY_ET_UNARY(name,functor)                            \
                                                                           \
template <typename T1>                                                     \
_bz_inline_et                                                              \
typename BZ_BLITZ_SCOPE(BzUnaryExprResult)<functor,T1>::T_result           \
name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1)                                 \
{									\
    typedef typename                                                       \
        BZ_BLITZ_SCOPE(BzUnaryExprResult)<functor,T1>::T_result result;    \
    return result(BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap()));       \
}

/*
 * Array expression templates: the macro BZ_DECLARE_ARRAY_ET_BINARY(X,Y)
 * declares a function or operator which takes two operands.
 * X is the function name (or operator), and Y is the functor object
 * which implements the operation.
 */

#define BZ_DECLARE_ARRAY_ET_BINARY(name, applic)                           \
                                                                           \
template <typename T1,typename T2>                                         \
_bz_inline_et                                                              \
typename BZ_BLITZ_SCOPE(BzBinaryExprResult)<applic,T1,T2>::T_result        \
name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1,                                 \
     const BZ_BLITZ_SCOPE(ETBase)<T2>& d2)                                 \
{                                                                          \
    typedef typename                                                       \
        BZ_BLITZ_SCOPE(BzBinaryExprResult)<applic,T1,T2>::T_result result; \
    return result(BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap()),        \
                  BZ_BLITZ_SCOPE(asExpr)<T2>::getExpr(d2.unwrap()));       \
}

#define BZ_DECLARE_ARRAY_ET_BINARY_TINYVEC(name, applic)                   \
                                                                           \
template <typename T1, typename T2, int N>                                 \
_bz_inline_et                                                              \
typename BZ_BLITZ_SCOPE(BzBinaryExprResult)<applic,                        \
    BZ_BLITZ_SCOPE(TinyVector)<T2,N>,T1>::T_result                         \
name(const BZ_BLITZ_SCOPE(TinyVector)<T2,N> d1,                            \
     const BZ_BLITZ_SCOPE(ETBase)<T1>& d2)                                 \
{                                                                          \
    typedef typename                                                       \
        BZ_BLITZ_SCOPE(BzBinaryExprResult)<applic,                         \
        BZ_BLITZ_SCOPE(TinyVector)<T2,N>,T1>::T_result result;             \
    return result(BZ_BLITZ_SCOPE(asExpr)<                                  \
                  BZ_BLITZ_SCOPE(TinyVector)<T2,N> >::getExpr(d1),         \
                  BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d2.unwrap()));       \
}                                                                          \
                                                                           \
template <typename T1, typename T2, int N>                                 \
_bz_inline_et                                                              \
typename BZ_BLITZ_SCOPE(BzBinaryExprResult)<applic,T1,                     \
    BZ_BLITZ_SCOPE(TinyVector)<T2,N> >::T_result                           \
name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1,                                 \
     const BZ_BLITZ_SCOPE(TinyVector)<T2,N> d2)                            \
{                                                                          \
    typedef typename                                                       \
        BZ_BLITZ_SCOPE(BzBinaryExprResult)<applic,T1,                      \
        BZ_BLITZ_SCOPE(TinyVector)<T2,N> >::T_result result;               \
    return result(BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap()),        \
                  BZ_BLITZ_SCOPE(asExpr)<                                  \
                  BZ_BLITZ_SCOPE(TinyVector)<T2,N> >::getExpr(d2));        \
}

#define BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(name, applic, sca)               \
                                                                           \
template<typename T>                                                       \
_bz_inline_et                                                              \
typename BZ_BLITZ_SCOPE(BzBinaryExprResult)<applic,sca,T>::T_result        \
name(const sca d1, const BZ_BLITZ_SCOPE(ETBase)<T>& d2)                    \
{                                                                          \
    typedef typename                                                       \
        BZ_BLITZ_SCOPE(BzBinaryExprResult)<applic,sca,T>::T_result result; \
    return result(BZ_BLITZ_SCOPE(asExpr)<sca >::getExpr(d1),               \
                  BZ_BLITZ_SCOPE(asExpr)<T>::getExpr(d2.unwrap()));        \
}                                                                          \
                                                                           \
template<typename T>                                                       \
_bz_inline_et                                                              \
typename BZ_BLITZ_SCOPE(BzBinaryExprResult)<applic,T,sca >::T_result       \
name(const BZ_BLITZ_SCOPE(ETBase)<T>& d1, const sca d2)                    \
{                                                                          \
    typedef typename                                                       \
        BZ_BLITZ_SCOPE(BzBinaryExprResult)<applic,T,sca >::T_result result;\
    return result(BZ_BLITZ_SCOPE(asExpr)<T>::getExpr(d1.unwrap()),         \
                  BZ_BLITZ_SCOPE(asExpr)<sca >::getExpr(d2));              \
}

/*
 * Array expression templates: the macro BZ_DECLARE_ARRAY_ET_TERNARY(X,Y)
 * declares a function or operator which takes three operands.
 * X is the function name (or operator), and Y is the functor object
 * which implements the operation.
 */

#define BZ_DECLARE_ARRAY_ET_TERNARY(name, applic)                          \
                                                                           \
template <typename T1, typename T2, typename T3>                           \
_bz_inline_et                                                              \
typename BZ_BLITZ_SCOPE(BzTernaryExprResult)<applic, T1, T2, T3>::T_result \
name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1,                                 \
     const BZ_BLITZ_SCOPE(ETBase)<T2>& d2,                                 \
     const BZ_BLITZ_SCOPE(ETBase)<T3>& d3)                                 \
{                                                                          \
    typedef typename                                                       \
        BZ_BLITZ_SCOPE(BzTernaryExprResult)<applic,T1,T2,T3>::T_result     \
        result;                                                            \
    return result(BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap()),        \
                  BZ_BLITZ_SCOPE(asExpr)<T2>::getExpr(d2.unwrap()),        \
                  BZ_BLITZ_SCOPE(asExpr)<T3>::getExpr(d3.unwrap()));       \
}


/*
 * Array expression templates: the macro BZ_DECLARE_ARRAY_ET_QUATERNARY(X,Y)
 * declares a function or operator which takes four operands.
 * X is the function name (or operator), and Y is the functor object
 * which implements the operation.
 */

#define BZ_DECLARE_ARRAY_ET_QUATERNARY(name, applic)			\
  									\
  template <typename T1, typename T2, typename T3, typename T4>		\
  _bz_inline_et								\
  typename BZ_BLITZ_SCOPE(BzQuaternaryExprResult)<applic, T1, T2, T3, T4>::T_result \
    name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1,				\
     const BZ_BLITZ_SCOPE(ETBase)<T2>& d2,				\
     const BZ_BLITZ_SCOPE(ETBase)<T3>& d3,				\
     const BZ_BLITZ_SCOPE(ETBase)<T4>& d4)				\
  {									\
  typedef typename							\
      BZ_BLITZ_SCOPE(BzQuaternaryExprResult)<applic,T1,T2,T3, T4>::T_result \
      result;								\
    return result(BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap()),	\
                  BZ_BLITZ_SCOPE(asExpr)<T2>::getExpr(d2.unwrap()),	\
                  BZ_BLITZ_SCOPE(asExpr)<T3>::getExpr(d3.unwrap()),	\
                  BZ_BLITZ_SCOPE(asExpr)<T4>::getExpr(d4.unwrap()));	\
}



#else /* !BZ_HAVE_TEMPLATES_AS_TEMPLATE_ARGUMENTS */

/*
 * Unary functions and operators
 */

#define BZ_DECLARE_ARRAY_ET_UNARY(name, functor)                           \
                                                                           \
template<typename T1>                                                      \
_bz_inline_et                                                              \
BZ_BLITZ_SCOPE(_bz_ArrayExpr)<BZ_BLITZ_SCOPE(_bz_ArrayExprUnaryOp)<        \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,                       \
    functor<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype> > >\
name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1)                                 \
{                                                                          \
    return BZ_BLITZ_SCOPE(_bz_ArrayExpr)<                                  \
        BZ_BLITZ_SCOPE(_bz_ArrayExprUnaryOp)<                              \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,                   \
        functor<_bz_typename                                               \
        BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype> > >(                \
        BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap()));                 \
}

/*
 * Array expression templates: the macro BZ_DECLARE_ARRAY_ET_BINARY(X,Y)
 * declares a function or operator which takes two operands.
 * X is the function name (or operator), and Y is the functor object
 * which implements the operation.
 */

#define BZ_DECLARE_ARRAY_ET_BINARY(name, applic)                           \
                                                                           \
template<typename T1, typename T2>                                         \
_bz_inline_et                                                              \
BZ_BLITZ_SCOPE(_bz_ArrayExpr)<BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)<       \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,                       \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr,                       \
    applic<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype,     \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_numtype> > >        \
name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1,                                 \
     const BZ_BLITZ_SCOPE(ETBase)<T2>& d2)                                 \
{                                                                          \
    return BZ_BLITZ_SCOPE(_bz_ArrayExpr)<                                  \
        BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)<                             \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,                   \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr,                   \
        applic<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype, \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_numtype> > >(   \
        BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap()),                  \
        BZ_BLITZ_SCOPE(asExpr)<T2>::getExpr(d2.unwrap()));                 \
}

#define BZ_DECLARE_ARRAY_ET_BINARY_TINYVEC(name, applic)                   \
                                                                           \
template <typename T1, typename T2, int N>                                 \
_bz_inline_et                                                              \
BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)<                                     \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<                                   \
    BZ_BLITZ_SCOPE(TinyVector)<T2,N> >::T_expr,                            \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,                       \
    applic<BZ_BLITZ_SCOPE(TinyVector)<T2,N>,                               \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype> >          \
name(const BZ_BLITZ_SCOPE(TinyVector)<T2,N> d1,                            \
     const BZ_BLITZ_SCOPE(ETBase)<T1>& d2)                                 \
{                                                                          \
    return BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)<                          \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<                               \
        BZ_BLITZ_SCOPE(TinyVector)<T2,N> >::T_expr,                        \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,                   \
        applic<BZ_BLITZ_SCOPE(TinyVector)<T2,N>,                           \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype> >(     \
        BZ_BLITZ_SCOPE(asExpr)<                                            \
        BZ_BLITZ_SCOPE(TinyVector)<T2,N> >::getExpr(d1),                   \
        BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d2.unwrap()));                 \
}                                                                          \
                                                                           \
template <typename T1, typename T2, int N>                                 \
_bz_inline_et                                                              \
BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)<                                     \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,                       \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<                                   \
    BZ_BLITZ_SCOPE(TinyVector)<T2,N> >::T_expr,                            \
    applic<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype,     \
    BZ_BLITZ_SCOPE(TinyVector)<T2,N> > >                                   \
name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1,                                 \
     const BZ_BLITZ_SCOPE(TinyVector)<T2,N> d2)                            \
{                                                                          \
    return BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)<                          \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,                   \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<                               \
        BZ_BLITZ_SCOPE(TinyVector)<T2,N> >::T_expr,                        \
        applic<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype, \
        BZ_BLITZ_SCOPE(TinyVector)<T2,N> > >(                              \
        BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap()),                  \
        BZ_BLITZ_SCOPE(asExpr)<                                            \
        BZ_BLITZ_SCOPE(TinyVector)<T2,N> >::getExpr(d2));                  \
}

#define BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(name, applic, sca)               \
                                                                           \
template<typename T>                                                       \
_bz_inline_et                                                              \
BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)<                                     \
    BZ_BLITZ_SCOPE(asExpr)<sca >::T_expr,                                  \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T>::T_expr,                        \
    applic<sca,_bz_typename BZ_BLITZ_SCOPE(asExpr)<T>::T_expr::T_numtype> >\
name(const sca d1, const BZ_BLITZ_SCOPE(ETBase)<T>& d2)                    \
{                                                                          \
    return BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)<                          \
        BZ_BLITZ_SCOPE(asExpr)<sca >::T_expr,                              \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T>::T_expr,                    \
        applic<sca,                                                        \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T>::T_expr::T_numtype> >(      \
        BZ_BLITZ_SCOPE(asExpr)<sca >::getExpr(d1),                         \
        BZ_BLITZ_SCOPE(asExpr)<T>::getExpr(d2.unwrap()));                  \
}                                                                          \
                                                                           \
template<typename T>                                                       \
_bz_inline_et                                                              \
BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)<                                     \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T>::T_expr,                        \
    BZ_BLITZ_SCOPE(asExpr)<sca >::T_expr,                                  \
    applic<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T>::T_expr::T_numtype,      \
           sca > >                                                         \
name(const BZ_BLITZ_SCOPE(ETBase)<T>& d1, const sca d2)                    \
{                                                                          \
    return BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)<                          \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T>::T_expr,                    \
        BZ_BLITZ_SCOPE(asExpr)<sca >::T_expr,                              \
        applic<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T>::T_expr::T_numtype,  \
               sca > >(                                                    \
        BZ_BLITZ_SCOPE(asExpr)<T>::getExpr(d1.unwrap()),                   \
        BZ_BLITZ_SCOPE(asExpr)<sca >::getExpr(d2));                        \
}

/*
 * Array expression templates: the macro BZ_DECLARE_ARRAY_ET_TERNARY(X,Y)
 * declares a function or operator which takes three operands.
 * X is the function name (or operator), and Y is the functor object
 * which implements the operation.
 */

#define BZ_DECLARE_ARRAY_ET_TERNARY(name, applic)                          \
                                                                           \
template<typename T1, typename T2, typename T3>                            \
_bz_inline_et                                                              \
BZ_BLITZ_SCOPE(_bz_ArrayExpr)<BZ_BLITZ_SCOPE(_bz_ArrayExprTernaryOp)<      \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,                       \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr,                       \
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T3>::T_expr,                       \
    applic<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype,     \
           _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_numtype,     \
           _bz_typename BZ_BLITZ_SCOPE(asExpr)<T3>::T_expr::T_numtype> > > \
name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1,                                 \
     const BZ_BLITZ_SCOPE(ETBase)<T2>& d2,                                 \
     const BZ_BLITZ_SCOPE(ETBase)<T3>& d3)                                 \
{                                                                          \
    return BZ_BLITZ_SCOPE(_bz_ArrayExpr)<                                  \
        BZ_BLITZ_SCOPE(_bz_ArrayExprTernaryOp)<                            \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,                   \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr,                   \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T3>::T_expr,                   \
        applic<_bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype, \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_numtype,        \
        _bz_typename BZ_BLITZ_SCOPE(asExpr)<T3>::T_expr::T_numtype> > >(   \
        BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap()),                  \
        BZ_BLITZ_SCOPE(asExpr)<T2>::getExpr(d2.unwrap()),                  \
        BZ_BLITZ_SCOPE(asExpr)<T3>::getExpr(d3.unwrap()));                 \
}


/*
 * Array expression templates: the macro BZ_DECLARE_ARRAY_ET_BINARY(X,Y)
 * declares a function or operator which takes two operands.
 * X is the function name (or operator), and Y is the functor object
 * which implements the operation.
 */

#define BZ_DECLARE_ARRAY_ET_QUATERNARY(name, functor)			\
  									\
  template<typename T1, typename T2, typename T3,typename T4>		\
  _bz_inline_et								\
  BZ_BLITZ_SCOPE(_bz_ArrayExpr)						\
    <									\
    BZ_BLITZ_SCOPE(_bz_ArrayExprQuaternaryOp)				\
    <									\
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,			\
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr,			\
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T3>::T_expr,			\
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T4>::T_expr,			\
    functor<								\
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype,		\
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_numtype,		\
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T3>::T_expr::T_numtype,		\
    _bz_typename BZ_BLITZ_SCOPE(asExpr)<T4>::T_expr::T_numtype		\
    > > >								\
									\
    name(const BZ_BLITZ_SCOPE(ETBase)<T1>& d1,				\
	 const BZ_BLITZ_SCOPE(ETBase)<T2>& d2,				\
	 const BZ_BLITZ_SCOPE(ETBase)<T3>& d3,				\
	 const BZ_BLITZ_SCOPE(ETBase)<T4>& d4)				\
  {									\
    return BZ_BLITZ_SCOPE(_bz_ArrayExpr)				\
      <									\
      BZ_BLITZ_SCOPE(_bz_ArrayExprBinaryOp)				\
      <									\
      _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr,			\
      _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr,			\
      functor								\
      <									\
      _bz_typename BZ_BLITZ_SCOPE(asExpr)<T1>::T_expr::T_numtype,	\
      _bz_typename BZ_BLITZ_SCOPE(asExpr)<T2>::T_expr::T_numtype,	\
      _bz_typename BZ_BLITZ_SCOPE(asExpr)<T3>::T_expr::T_numtype,	\
      _bz_typename BZ_BLITZ_SCOPE(asExpr)<T4>::T_expr::T_numtype	\
      > > >								\
      (									\
       BZ_BLITZ_SCOPE(asExpr)<T1>::getExpr(d1.unwrap()),		\
       BZ_BLITZ_SCOPE(asExpr)<T2>::getExpr(d2.unwrap()),		\
       BZ_BLITZ_SCOPE(asExpr)<T2>::getExpr(d3.unwrap()),		\
       BZ_BLITZ_SCOPE(asExpr)<T2>::getExpr(d4.unwrap()));		\
}


#endif /* BZ_HAVE_TEMPLATES_AS_TEMPLATE_ARGUMENTS */

/*
 * User-defined expression template routines
 */

#define BZ_DECLARE_FUNCTION(name)                                          \
BZ_DEFINE_UNARY_FUNC(name ## _impl,name,true)				\
BZ_DECLARE_ARRAY_ET_UNARY(name,name ## _impl)

#define BZ_DECLARE_FUNCTION_RET(name,return_type)                          \
BZ_DEFINE_UNARY_FUNC_RET(name ## _impl,name,return_type)                   \
BZ_DECLARE_ARRAY_ET_UNARY(name,name ## _impl)

#define BZ_DECLARE_FUNCTION2(name)                                         \
BZ_DEFINE_BINARY_FUNC(name ## _impl,name)                                  \
BZ_DECLARE_ARRAY_ET_BINARY(name, name ## _impl)

#define BZ_DECLARE_FUNCTION2_RET(name,return_type)                         \
BZ_DEFINE_BINARY_FUNC_RET(name ## _impl,name,return_type)                  \
BZ_DECLARE_ARRAY_ET_BINARY(name, name ## _impl)

#define BZ_DECLARE_FUNCTION2_SCALAR(name, sca)                             \
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(name, name ## _impl, sca)

#define BZ_DECLARE_FUNCTION3(name)                                         \
BZ_DEFINE_TERNARY_FUNC(name ## _impl,name)                                 \
BZ_DECLARE_ARRAY_ET_TERNARY(name, name ## _impl)

#define BZ_DECLARE_FUNCTION3_RET(name,return_type)                         \
BZ_DEFINE_TERNARY_FUNC_RET(name ## _impl,name,return_type)                 \
BZ_DECLARE_ARRAY_ET_TERNARY(name, name ## _impl)

#define BZ_DECLARE_FUNCTION4(name)					\
  BZ_DEFINE_QUATERNARY_FUNC(name ## _impl,name)				\
  BZ_DECLARE_ARRAY_ET_QUATERNARY(name, name ## _impl)


BZ_NAMESPACE_END

#endif








