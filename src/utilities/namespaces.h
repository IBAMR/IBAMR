#ifndef included_namespaces
#define included_namespaces

// Filename: namespaces.h
// Last modified: <27.Jun.2010 15:34:36 griffith@griffith-macbook-pro.local>
// Created on 27 Jun 2010 by Boyce Griffith (griffith@boyce-griffiths-macbook-pro.local)

//////////////////////////////////////////////////////////////////////////////

/*!
 * Defines "using" declarations for all SAMRAI and IBTK namespaces.  This header
 * file may be included in application codes, but it MUST NOT be included in any
 * other header (.h) or inline (.I) file in the library.
 */
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

namespace IBTK
{
}
using namespace IBTK;

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_namespaces
