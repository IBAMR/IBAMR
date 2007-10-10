// Filename: IBMarkerRefineOperator.C
// Last modified: <10.Oct.2007 00:44:50 griffith@box221.cims.nyu.edu>
// Created on 04 Oct 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBMarkerRefineOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/IBMarker.h>

// STOOLS INCLUDES
#include <stools/STOOLS_Utilities.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <IndexData.h>
#include <IndexVariable.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string IBMarkerRefineOperator::s_op_name = "IB_MARKER_REFINE";

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBMarkerRefineOperator::IBMarkerRefineOperator()
{
    // intentionally blank
    return;
}// IBMarkerRefineOperator

IBMarkerRefineOperator::~IBMarkerRefineOperator()
{
    // intentionally blank
    return;
}// ~IBMarkerRefineOperator

bool
IBMarkerRefineOperator::findRefineOperator(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& var,
    const std::string& op_name) const
{
    SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<NDIM,IBMarker> > mark_var = var;
    return (!mark_var.isNull() && op_name == s_op_name);
}// findRefineOperator

const std::string&
IBMarkerRefineOperator::getOperatorName() const
{
    return s_op_name;
}// getOperatorName

int
IBMarkerRefineOperator::getOperatorPriority() const
{
    return 0;
}// getOperatorPriority

SAMRAI::hier::IntVector<NDIM>
IBMarkerRefineOperator::getStencilWidth() const
{
    return 2;
}// getStencilWidth

void
IBMarkerRefineOperator::refine(
    SAMRAI::hier::Patch<NDIM>& fine,
    const SAMRAI::hier::Patch<NDIM>& coarse,
    const int dst_component,
    const int src_component,
    const SAMRAI::hier::Box<NDIM>& fine_box,
    const SAMRAI::hier::IntVector<NDIM>& ratio) const
{
    SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBMarker> > dst_mark_data = fine  .getPatchData(dst_component);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBMarker> > src_mark_data = coarse.getPatchData(src_component);

    const SAMRAI::hier::Box<NDIM>& fine_patch_box = fine.getBox();
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > fine_patch_geom =
        fine.getPatchGeometry();
    const SAMRAI::hier::Index<NDIM>& fine_patch_lower = fine_patch_box.lower();
    const SAMRAI::hier::Index<NDIM>& fine_patch_upper = fine_patch_box.upper();
    const double* const fine_patchXLower = fine_patch_geom->getXLower();
    const double* const fine_patchXUpper = fine_patch_geom->getXUpper();
    const double* const fine_patchDx = fine_patch_geom->getDx();

    const SAMRAI::hier::Box<NDIM> coarse_box = SAMRAI::hier::Box<NDIM>::coarsen(fine_box,ratio);
    for (SAMRAI::pdat::IndexData<NDIM,IBMarker>::Iterator it(*src_mark_data); it; it++)
    {
        const SAMRAI::hier::Index<NDIM>& coarse_i = it.getIndex();
        if (coarse_box.contains(coarse_i))
        {
            const IBMarker& coarse_mark = it();
            const std::vector<double>& coarse_X = coarse_mark.getPositions();
            const std::vector<double>& coarse_U = coarse_mark.getVelocities();
            const std::vector<int>& coarse_idx = coarse_mark.getIndices();
            for (int k = 0; k < coarse_mark.getNumberOfMarkers(); ++k)
            {
                const double* const X = &coarse_X[NDIM*k];
                const double* const U = &coarse_U[NDIM*k];
                const int& idx = coarse_idx[k];
                const SAMRAI::hier::Index<NDIM> fine_i =
                    STOOLS::STOOLS_Utilities::getCellIndex(
                        X,fine_patchXLower,fine_patchXUpper,fine_patchDx,fine_patch_lower,fine_patch_upper);
                if (fine_box.contains(fine_i))
                {
                    if (!dst_mark_data->isElement(fine_i))
                    {
                        dst_mark_data->appendItem(fine_i, IBMarker());
                    }
                    IBMarker& fine_mark = *(dst_mark_data->getItem(fine_i));
                    std::vector<double>& fine_X = fine_mark.getPositions();
                    std::vector<double>& fine_U = fine_mark.getVelocities();
                    std::vector<int>& fine_idx = fine_mark.getIndices();

                    fine_X.insert(fine_X.end(),X,X+NDIM);
                    fine_U.insert(fine_U.end(),U,U+NDIM);
                    fine_idx.push_back(idx);
                }
            }
        }
    }
    return;
}// refine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBMarkerRefineOperator>;

//////////////////////////////////////////////////////////////////////////////
