// Filename: INSStaggeredPhysicalBoundaryHelper.C
// Last modified: <22.Jul.2008 19:10:26 griffith@box230.cims.nyu.edu>
// Created on 22 Jul 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredPhysicalBoundaryHelper.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <SideIndex.h>
#include <SideData.h>
#include <tbox/MathUtilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredPhysicalBoundaryHelper::INSStaggeredPhysicalBoundaryHelper()
    : d_hierarchy(NULL),
      d_dirichlet_bdry()
{
    // intentionally blank
    return;
}// INSStaggeredPhysicalBoundaryHelper

INSStaggeredPhysicalBoundaryHelper::~INSStaggeredPhysicalBoundaryHelper()
{
    // intentionally blank
    return;
}// ~INSStaggeredPhysicalBoundaryHelper

void
INSStaggeredPhysicalBoundaryHelper::zeroValuesAtDirichletBoundaries(
    const int patch_data_idx,
    const int coarsest_ln,
    const int finest_ln) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln);
         ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(patch_num);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > patch_data = patch->getPatchData(patch_data_idx);

            // Compute the boundary fill boxes.
            const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> > physical_codim1_boxes = IBTK::PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
            const int n_physical_codim1_boxes = physical_codim1_boxes.size();

            // Compute the locations of the Dirichlet boundary.
            const std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,int> > >& dirichlet_bdry = (*d_dirichlet_bdry[ln].find(patch_num)).second;
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const int location_index   = bdry_box.getLocationIndex();
                const int bdry_normal_axis = location_index / 2;
                const SAMRAI::hier::Box<NDIM>& bc_coef_box = dirichlet_bdry[n]->getBox();
                for (SAMRAI::hier::Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = it();
                    if ((*dirichlet_bdry[n])(i,0) == 1)
                    {
                        const SAMRAI::pdat::SideIndex<NDIM> i_s(i, bdry_normal_axis, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                        (*patch_data)(i_s) = 0.0;
                    }
                }
            }
        }
    }
    return;
}// zeroValuesAtDirichletBoundaries

void
INSStaggeredPhysicalBoundaryHelper::cacheBcCoefData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& u_var,
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const double fill_time,
    const SAMRAI::hier::IntVector<NDIM>& gcw_to_fill,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >& hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
#endif
    if (!d_hierarchy.isNull()) clearBcCoefData();
    d_hierarchy = hierarchy;
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    d_dirichlet_bdry.resize(finest_hier_level+1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(patch_num);
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

            // Compute the boundary fill boxes.
            const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> > physical_codim1_boxes = IBTK::PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
            const int n_physical_codim1_boxes = physical_codim1_boxes.size();

            // Compute the locations of the Dirichlet boundary.
            std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,int> > >& dirichlet_bdry = d_dirichlet_bdry[ln][patch_num];
            dirichlet_bdry.resize(n_physical_codim1_boxes,SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,int> >(NULL));
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const SAMRAI::hier::Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                const SAMRAI::hier::BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
                const int location_index   = bdry_box.getLocationIndex();
                const int bdry_normal_axis = location_index / 2;

                const SAMRAI::hier::Box<NDIM> bc_coef_box = IBTK::PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > acoef_data = new SAMRAI::pdat::ArrayData<NDIM,double>(bc_coef_box, 1);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > bcoef_data = new SAMRAI::pdat::ArrayData<NDIM,double>(bc_coef_box, 1);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > gcoef_data = NULL;
                u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, u_var, *patch, trimmed_bdry_box, fill_time);

                dirichlet_bdry[n] = new SAMRAI::pdat::ArrayData<NDIM,int>(bc_coef_box, 1);
                for (SAMRAI::hier::Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = it();
                    const double& alpha = (*acoef_data)(i,0);
                    const double& beta  = (*bcoef_data)(i,0);
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(SAMRAI::tbox::MathUtilities<double>::equalEps(alpha+beta,1.0));
                    TBOX_ASSERT(SAMRAI::tbox::MathUtilities<double>::equalEps(alpha,1.0) || SAMRAI::tbox::MathUtilities<double>::equalEps(beta,1.0));
#endif
                    (*dirichlet_bdry[n])(i,0) = SAMRAI::tbox::MathUtilities<double>::equalEps(alpha,1.0) || !SAMRAI::tbox::MathUtilities<double>::equalEps(beta,1.0);
                }
            }
        }
    }
    return;
}// cacheBcCoefData

void
INSStaggeredPhysicalBoundaryHelper::clearBcCoefData()
{
    d_hierarchy.setNull();
    d_dirichlet_bdry.clear();
    return;
}// clearBcCoefData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredPhysicalBoundaryHelper>;

//////////////////////////////////////////////////////////////////////////////
