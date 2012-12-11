// -*- C++ -*-
/***************************************************************************
 * blitz/bzconfig.h      Select compiler-specific config file
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

/* Select the compiler-specific config.h header file */
#ifndef BZCONFIG_H
#define BZCONFIG_H

#if defined(__APPLE__)
#if defined(__GNUC__)
/* GNU gcc compiler for newer Mac OS X Darwin */
#include <blitz/gnu/bzconfig.h>
#else
/* IBM xlc compiler for Darwin */
#include <blitz/apple/bzconfig.h>
#endif

#elif defined(__INTEL_COMPILER)
/* Intel icc compiler */
#include <blitz/intel/bzconfig.h>

#elif defined(_MSC_VER)
/* Microsoft VS.NET compiler */
#include <blitz/ms/bzconfig.h>

#elif defined(__xlC__)
/* IBM xlC compiler */
#include <blitz/ibm/bzconfig.h>

#elif defined(__DECCXX)
/* Compaq cxx compiler */
#include <blitz/compaq/bzconfig.h>

#elif defined(__HP_aCC)
/* HP aCC compiler */
#include <blitz/hp/bzconfig.h>

#elif defined(_SGI_COMPILER_VERSION)
/* SGI CC compiler */
#include <blitz/sgi/bzconfig.h>

#elif defined(__SUNPRO_CC)
/* SunPRO CC compiler */
#include <blitz/sun/bzconfig.h>

#elif defined(__PATHCC__)
/* Pathscale pathCC compiler */
#include <blitz/pathscale/bzconfig.h>

#elif defined(__GNUC__)
/* GNU gcc compiler */
#include <blitz/gnu/bzconfig.h>

#elif defined(__PGI)
/* PGI pgCC compiler */
#include <blitz/pgi/bzconfig.h>

#elif defined(__KCC)
/* KAI KCC compiler */
#include <blitz/kai/bzconfig.h>

#elif defined(__FUJITSU)
/* Fujitsu FCC compiler */
#include <blitz/fujitsu/bzconfig.h>

/* Add other compilers here */

#else
#error Unknown compiler
#endif

#endif /* BZCONFIG_H */
