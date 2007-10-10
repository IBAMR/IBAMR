// Filename: IBMarkerCoarsenOperator.C
// Last modified: <10.Oct.2007 00:10:39 griffith@box221.cims.nyu.edu>
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

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <IndexData.h>
#include <IndexVariable.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string IBMarkerCoarsenOperator::s_op_name = "IB_MARKER_COARSEN";

namespace
{
inline int
coarsen(
    const int index,
    const int ratio)
{
    return (index < 0 ? (index+1)/ratio-1 : index/ratio);
}// coarsen

inline SAMRAI::hier::Index<NDIM>
coarsen_index(
    const SAMRAI::hier::Index<NDIM>& i,
    const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    SAMRAI::hier::Index<NDIM> coarse_i;
    for (int d = 0; d < NDIM; ++d)
    {
        coarse_i(d) = coarsen(i(d), ratio(d));
    }
    return coarse_i;
}// coarsen_index
}

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

    const SAMRAI::hier::Box<NDIM> fine_box = SAMRAI::hier::Box<NDIM>::refine(coarse_box,ratio);
    for (SAMRAI::pdat::IndexData<NDIM,IBMarker>::Iterator it(*src_mark_data); it; it++)
    {
        const SAMRAI::hier::Index<NDIM>& fine_i = it.getIndex();
        const SAMRAI::hier::Index<NDIM> coarse_i = coarsen_index(fine_i,ratio);
        if (fine_box.contains(fine_i) && coarse_box.contains(coarse_i))
        {
            const IBMarker& fine_mark = it();
            const std::vector<double>& fine_X = fine_mark.getPositions();
            const std::vector<double>& fine_U = fine_mark.getVelocities();
            const std::vector<int>& fine_idx = fine_mark.getIndices();

            if (!dst_mark_data->isElement(coarse_i))
            {
                dst_mark_data->appendItem(coarse_i, IBMarker());
            }
            IBMarker& coarse_mark = *(dst_mark_data->getItem(coarse_i));
            std::vector<double>& coarse_X = coarse_mark.getPositions();
            std::vector<double>& coarse_U = coarse_mark.getVelocities();
            std::vector<int>& coarse_idx = coarse_mark.getIndices();

            coarse_X.insert(coarse_X.end(),fine_X.begin(),fine_X.end());
            coarse_U.insert(coarse_U.end(),fine_U.begin(),fine_U.end());
            coarse_idx.insert(coarse_idx.end(),fine_idx.begin(),fine_idx.end());
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
