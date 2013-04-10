/***************************************************************************
 * blitz/tinyvec.cc  Declaration of TinyVector methods
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
 ***************************************************************************/

#ifndef BZ_TVEVALUATE_H
#define BZ_TVEVALUATE_H

#include <blitz/tinyvec2.h>
#include <blitz/update.h>
#include <blitz/blitz.h>
#include <blitz/meta/vecassign.h>

BZ_NAMESPACE(blitz)


/** The _tv_evaluator class has a bool template argument that is used
    to select code paths at compile time.  */
template<bool unroll, int N_length>
struct _tv_evaluator {

  /** The select_evaluation function redirects expressions that do not
      contains solely TinyVector operands to the general evaluation
      function. The generic template (for unroll=false, note that
      "unroll" us the wrong name for this function, the template
      parameter in this context really means "use_full_eval") uses the
      TinyVector-only evaluation. Since TinyVectors can't have funny
      storage, ordering, stride, or anything, it's now just a matter
      of evaluating it like in the old vecassign. */
  template<typename T, typename T_expr, typename T_update>
  static _bz_forceinline void
  select_evaluation(TinyVector<T, N_length>& dest,
		    const T_expr& expr, T_update) {

    // since we can't resize tinyvectors, there are two options: all
    // vectors have our size or the expression is malformed.
    // Check that all operands have the same shape
#ifdef BZ_DEBUG
    if (!expr.shapeCheck(dest.shape()))
      {
	if (assertFailMode == false)
	  {
	    cerr << "[Blitz++] Shape check failed: Module " << __FILE__
		 << " line " << __LINE__ << endl
		 << "          Expression: ";
        prettyPrintFormat format(true);   // Use terse formatting
        BZ_STD_SCOPE(string) str;
        expr.prettyPrint(str, format);
        cerr << str << endl ;
      }
    }
#endif

  BZPRECHECK(expr.shapeCheck(dest.shape()),
	     "Shape check failed." << endl << "Expression:");

  BZPRECONDITION(expr.isUnitStride(0));
  BZPRECONDITION(T_expr::rank_<=1);
  BZPRECONDITION(T_expr::numIndexPlaceholders==0);

  // now call the aligned (unrolled or not) evaluation function
  const bool do_unroll = N_length < BZ_TV_EVALUATE_UNROLL_LENGTH;
  _tv_evaluator<do_unroll, N_length>::evaluate_aligned(dest.data(), expr, T_update());
  }

  /** This version of the evaluation function assumes that the
      TinyVectors have appropriate alignment (as will always be the
      case if they are actual TinyVector objects and not created using
      reinterpret_cast in the chunked_updater. If no SIMD width is
      set, however, we can not guarantee alignment and must remove the
      pragmas. */
  template<typename T_numtype, typename T_expr, typename T_update>
  static _bz_forceinline void
  evaluate_aligned(T_numtype* data, const T_expr& expr, T_update) {
#ifdef BZ_USE_ALIGNMENT_PRAGMAS
#pragma ivdep
#pragma vector aligned
#endif
    for (int i=0; i < N_length; ++i)
      T_update::update(data[i], expr.fastRead(i));
  }

  /** This version of the evaluation function is used when vectorizing
      expressions that we know can't be aligned. The only difference
      with evaluate_aligned is the compiler pragma that tells the
      compiler it is unaligned. */
  template<typename T_numtype, typename T_expr, typename T_update>
  static _bz_forceinline void
  evaluate_unaligned(T_numtype* data, const T_expr& expr, T_update) {
#ifdef BZ_USE_ALIGNMENT_PRAGMAS
#pragma ivdep
#pragma vector unaligned
#endif
  for (int i=0; i < N_length; ++i)
    T_update::update(data[i], expr.fastRead(i));
  }
};

/** Specialization of the _tv_evaluator class for false template arguments. */
template<int N_length>
struct _tv_evaluator<true, N_length> {

  /** The false version of select_evaluation is picked for expressions
      that contain operands other than TinyVectors. It just redirects
      to the general evaluation function. */
  template<typename T, typename T_expr, typename T_update>
  static _bz_forceinline void
  select_evaluation(TinyVector<T, N_length>& dest,
		    const T_expr& expr, T_update) {
    _bz_evaluate(dest, expr, T_update());
  }

  /** This version of the evaluation function assumes that the
      TinyVectors have appropriate alignment (as will always be the
      case if they are actual TinyVector objects and not created using
      reinterpret_cast in the chunked_updater. */
  template<typename T_numtype, typename T_expr, typename T_update>
  static _bz_forceinline void
  evaluate_aligned(T_numtype* data, const T_expr& expr, T_update) {
#ifdef BZ_USE_ALIGNMENT_PRAGMAS
    //#pragma ivdep
    //#pragma vector aligned
#endif
  _bz_meta_vecAssign<N_length, 0>::fastAssign(data, expr, T_update());
  }

  /** This version of the evaluation function is used when vectorizing
      expressions that we know can't be aligned. The only difference
      with evaluate_aligned is the compiler pragma that tells the
      compiler it is unaligned. */
  template<typename T_numtype, typename T_expr, typename T_update>
  static _bz_forceinline void
  evaluate_unaligned(T_numtype* data, const T_expr& expr, T_update) {
    //#pragma ivdep
    //#pragma vector unaligned
  _bz_meta_vecAssign<N_length, 0>::fastAssign(data, expr, T_update());
  }
};


/** This function selects evaluation path by calling select_evaluation
    with a bool argument which is false if the expression only
    contains TinyVector operands. */
template<typename P_numtype, int N_length>
template<typename T_expr, typename T_update>
_bz_forceinline
void
TinyVector<P_numtype,N_length>::_tv_evaluate(const T_expr& expr, T_update)
{
  const bool mixed_expr =
    (T_expr::numArrayOperands>0) ||
    (T_expr::numTMOperands>0) ||
    (T_expr::numIndexPlaceholders>0);
  _tv_evaluator<mixed_expr, N_length>::select_evaluation(*this, expr, T_update());
}


BZ_NAMESPACE_END

#endif // BZ_TVEVALUATE_H
