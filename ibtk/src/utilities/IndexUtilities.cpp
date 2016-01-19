// Filename: IndexUtilities.cpp
// Created on 12 Jul 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "CartesianGridGeometry.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

Pointer<CartesianGridGeometry<NDIM> > IndexUtilities::s_grid_geometry = Pointer<CartesianGridGeometry<NDIM> >(NULL);
const double* IndexUtilities::s_dx0 = NULL;
const double* IndexUtilities::s_x_lower = NULL;
const double* IndexUtilities::s_x_upper = NULL;
IntVector<NDIM> IndexUtilities::s_ilower;
IntVector<NDIM> IndexUtilities::s_iupper;
bool IndexUtilities::s_is_initialized = false;

/////////////////////////////// PUBLIC ///////////////////////////////////////
void IndexUtilities::init(Pointer<CartesianGridGeometry<NDIM> > grid_geometry)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geometry);
#endif

    s_grid_geometry = grid_geometry;
    s_dx0 = s_grid_geometry->getDx();
    s_x_lower = s_grid_geometry->getXLower();
    s_x_upper = s_grid_geometry->getXUpper();

    const BoxArray<NDIM>& phys_domain = s_grid_geometry->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(phys_domain.size() == 1);
#endif
    s_ilower = phys_domain[0].lower();
    s_iupper = phys_domain[0].upper();

    s_is_initialized = true;

    return;
} // init

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
