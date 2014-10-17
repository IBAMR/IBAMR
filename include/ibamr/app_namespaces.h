// Filename: app_namespaces.h
// Created on 19 Aug 2011 by Boyce Griffith
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

#ifndef included_app_namespaces
#define included_app_namespaces

//////////////////////////////////////////////////////////////////////////////

/*!
 * Defines "using" declarations for all namespaces used in IBAMR and IBTK.  This
 * header file may be included in application codes, but it MUST NOT be included
 * in any other header (.h) or inline (.I) file in the library.
 */
namespace Eigen
{
}
using namespace Eigen;

namespace IBAMR
{
}
using namespace IBAMR;

namespace IBTK
{
}
using namespace IBTK;

namespace SAMRAI
{
namespace algs
{
}
namespace appu
{
}
namespace geom
{
}
namespace hier
{
}
namespace math
{
}
namespace mesh
{
}
namespace pdat
{
}
namespace solv
{
}
namespace tbox
{
}
namespace xfer
{
}
}
using namespace SAMRAI;
using namespace SAMRAI::algs;
using namespace SAMRAI::appu;
using namespace SAMRAI::geom;
using namespace SAMRAI::hier;
using namespace SAMRAI::math;
using namespace SAMRAI::mesh;
using namespace SAMRAI::pdat;
using namespace SAMRAI::solv;
using namespace SAMRAI::tbox;
using namespace SAMRAI::xfer;

namespace libMesh
{
}
using namespace libMesh;

namespace std
{
}
using namespace std;

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_app_namespaces
