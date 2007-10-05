// Filename: IBMarkerCoarsenOperator.C
// Last modified: <04.Oct.2007 23:49:51 griffith@box221.cims.nyu.edu>
// Created on 30 Sep 2006 by Boyce Griffith (boyce@trasnaform2.local)

#include "IBMarkerCoarsenOperator.h"

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

const std::string IBMarkerCoarsenOperator::s_op_name = "IB_MARKER_COARSEN";

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBMarkerCoarsenOperator::IBMarkerCoarsenOperator()
{
    // intentionally blank
    return;
}// IBMarkerCoarsenOperator

IBMarkerCoarsenOperator::~IBMarkerCoarsenOperator()
{
    // intentionally blank
    return;
}// ~IBMarkerCoarsenOperator

bool
IBMarkerCoarsenOperator::findCoarsenOperator(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& var,
    const std::string &op_name) const
{
    SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<NDIM,IBMarker> > mark_var = var;
    return (!mark_var.isNull() && op_name == s_op_name);
}// findCoarsenOperator

const std::string&
IBMarkerCoarsenOperator::getOperatorName() const
{
    return s_op_name;
}// getOperatorName

int
IBMarkerCoarsenOperator::getOperatorPriority() const
{
    return 128;
}// getOperatorPriority

SAMRAI::hier::IntVector<NDIM>
IBMarkerCoarsenOperator::getStencilWidth() const
{
    return 0;
}// getStencilWidth

void
IBMarkerCoarsenOperator::coarsen(
    SAMRAI::hier::Patch<NDIM>& coarse,
    const SAMRAI::hier::Patch<NDIM>& fine,
    const int dst_component,
    const int src_component,
    const SAMRAI::hier::Box<NDIM>& coarse_box,
    const SAMRAI::hier::IntVector<NDIM>& ratio) const
{
    SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBMarker> > dst_mark_data = coarse.getPatchData(dst_component);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBMarker> > src_mark_data = fine  .getPatchData(src_component);

    const SAMRAI::hier::Box<NDIM>& coarse_patch_box = coarse.getBox();
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > coarse_patch_geom =
        coarse.getPatchGeometry();
    const SAMRAI::hier::Index<NDIM>& coarse_patch_lower = coarse_patch_box.lower();
    const SAMRAI::hier::Index<NDIM>& coarse_patch_upper = coarse_patch_box.upper();
    const double* const coarse_patchXLower = coarse_patch_geom->getXLower();
    const double* const coarse_patchXUpper = coarse_patch_geom->getXUpper();
    const double* const coarse_patchDx = coarse_patch_geom->getDx();

    const SAMRAI::hier::Box<NDIM> fine_box = SAMRAI::hier::Box<NDIM>::refine(coarse_box,ratio);
    for (SAMRAI::pdat::IndexData<NDIM,IBMarker>::Iterator it(*src_mark_data); it; it++)
    {
        const SAMRAI::hier::Index<NDIM>& fine_i = it.getIndex();
        if (fine_box.contains(fine_i))
        {
            const IBMarker& fine_mark = it();

            const int num_marks = fine_mark.getNumberOfMarkers();
            const std::vector<double>& fine_X = fine_mark.getPositions();
            const std::vector<double>& fine_U = fine_mark.getVelocities();
            const std::vector<int>& fine_idx = fine_mark.getIndices();

            for (int k = 0; k < num_marks; ++k)
            {
                const double* const X = &fine_X[NDIM*k];
                const double* const U = &fine_U[NDIM*k];
                const int& idx = fine_idx[k];
                const SAMRAI::hier::Index<NDIM> coarse_i =
                    STOOLS::STOOLS_Utilities::getCellIndex(
                        X,coarse_patchXLower,coarse_patchXUpper,coarse_patchDx,coarse_patch_lower,coarse_patch_upper);
                if (coarse_box.contains(coarse_i))
                {
                    if (!dst_mark_data->isElement(coarse_i))
                    {
                        dst_mark_data->appendItem(coarse_i, IBMarker());
                    }
                    IBMarker* const coarse_mark = dst_mark_data->getItem(coarse_i);
                    coarse_mark->getPositions() .insert(coarse_mark->getPositions() .end(),X,X+NDIM);
                    coarse_mark->getVelocities().insert(coarse_mark->getVelocities().end(),U,U+NDIM);
                    coarse_mark->getIndices().push_back(idx);
                }
            }
        }
    }
    return;
}// coarsen

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBMarkerCoarsenOperator>;

//////////////////////////////////////////////////////////////////////////////
