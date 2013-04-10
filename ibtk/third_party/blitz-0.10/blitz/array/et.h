// -*- C++ -*-
/***************************************************************************
 * blitz/array/et.h  Include expression templates implementation for arrays
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
#ifndef BZ_ARRAY_ET_H
#define BZ_ARRAY_ET_H

#ifdef BZ_NEW_EXPRESSION_TEMPLATES
 #include <blitz/array/newet.h>     // Expression templates
#else
 #include <blitz/array/bops.cc>     // Expression templates, two operands
 #include <blitz/array/uops.cc>     // Expression templates, math functions
 #include <blitz/array/misc.cc>     // Expression templates, miscellaneous
#endif

#endif
