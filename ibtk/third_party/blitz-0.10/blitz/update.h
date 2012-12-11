// -*- C++ -*-
/***************************************************************************
 * blitz/update.h      Declaration of the _bz_XXXX updater classes
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

#ifndef BZ_UPDATE_H
#define BZ_UPDATE_H

#include <blitz/blitz.h>

BZ_NAMESPACE(blitz)

class _bz_updater_base { };

#define BZ_DECL_UPDATER(name,op,symbol)					\
  template<typename T_dest, typename T_source>                          \
  class name : public _bz_updater_base {				\
  public:								\
  template<typename T1, typename T2> struct updateCast {		\
    typedef name<T1, T2> T_updater;					\
  };									\
  									\
  static inline void update(T_dest& restrict x, T_source y)		\
    { x op y; }								\
  static void prettyPrint(BZ_STD_SCOPE(string) &str)			\
    { str += symbol; }							\
  }

template<typename T_dest, typename T_source>
class _bz_update : public _bz_updater_base {
  public:

  /** Traits class used to get an updater with different types. */
  template<typename T1, typename T2> struct updateCast {
    typedef _bz_update<T1, T2> T_updater;
  };

    static inline void update(T_dest& restrict x, T_source y)
  { x = /*(X)*/y; }

    static void prettyPrint(BZ_STD_SCOPE(string) &str)
    { str += "="; }
};

BZ_DECL_UPDATER(_bz_plus_update, +=, "+=");
BZ_DECL_UPDATER(_bz_minus_update, -=, "-=");
BZ_DECL_UPDATER(_bz_multiply_update, *=, "*=");
BZ_DECL_UPDATER(_bz_divide_update, /=, "/=");
BZ_DECL_UPDATER(_bz_mod_update, %=, "%=");
BZ_DECL_UPDATER(_bz_xor_update, ^=, "^=");
BZ_DECL_UPDATER(_bz_bitand_update, &=, "&=");
BZ_DECL_UPDATER(_bz_bitor_update, |=, "|=");
BZ_DECL_UPDATER(_bz_shiftl_update, <<=, "<<=");
BZ_DECL_UPDATER(_bz_shiftr_update, >>=, ">>=");

BZ_NAMESPACE_END

#endif // BZ_UPDATE_H

