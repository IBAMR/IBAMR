// -*- C++ -*-
/***************************************************************************
 * blitz/tau.h       Integration with the Tau profiling package.
 *                   See http://www.acl.lanl.gov/tau/
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

#ifndef BZ_TAU_H
#define BZ_TAU_H

#ifdef BZ_TAU_PROFILING
 #define TAU_BLITZ  TAU_USER1
 #include <Profile/Profiler.h>

#else
 #define TYPE_STRING(profileString, str)
 #define PROFILED_BLOCK(name, type)
 #define TAU_TYPE_STRING(profileString, str)
 #define TAU_PROFILE(name, type, group)
 #define TAU_PROFILE_TIMER(var, name, type, group)
 #define TAU_PROFILE_START(var)
 #define TAU_PROFILE_STOP(var)
 #define TAU_PROFILE_STMT(stmt)
 #define TAU_PROFILE_EXIT(msg)
 #define TAU_PROFILE_INIT(argc, argv)
 #define TAU_PROFILE_SET_NODE(node)
 #define CT(obj)
#endif // ! BZ_TAU_PROFILING

#endif // BZ_TAU_H
