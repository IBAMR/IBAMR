// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_CCScaledLaplaceOperator
#define included_IBTK_CCScaledLaplaceOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <SAMRAIVectorReal.h>

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Cell-centered operator A*u = C*u + S*div(D*grad(u)).
 *
 * C and D are specified through the inherited PoissonSpecifications.
 * S is stored in scalar cell-centered patch data.
 * The same S is used for all vector components and data depths.
 *
 * Ghost filling and boundary-condition handling are inherited from
 * CCLaplaceOperator.
 */
class CCScaledLaplaceOperator : public CCLaplaceOperator
{
public:
    CCScaledLaplaceOperator(std::string object_name,
                            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr,
                            bool homogeneous_bc = true);

    ~CCScaledLaplaceOperator() override = default;

    /*!
     * \brief Use scalar cell-centered patch data as the scale.
     *
     * The caller manages allocation and coefficient values on all levels.
     */
    void setScalePatchDataIndex(int scale_idx);

protected:
    void applySpecialized(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                          SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

private:
    CCScaledLaplaceOperator() = delete;
    CCScaledLaplaceOperator(const CCScaledLaplaceOperator&) = delete;
    CCScaledLaplaceOperator& operator=(const CCScaledLaplaceOperator&) = delete;

    int d_scale_idx = IBTK::invalid_index;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CCScaledLaplaceOperator
