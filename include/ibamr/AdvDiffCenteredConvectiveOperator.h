// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_AdvDiffCenteredConvectiveOperator
#define included_IBAMR_AdvDiffCenteredConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/CellConvectiveOperator.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffCenteredConvectiveOperator is a concrete
 * ConvectiveOperator that implements a centered convective differencing
 * operator.
 *
 * \see AdvDiffSemiImplicitHierarchyIntegrator
 */
class AdvDiffCenteredConvectiveOperator : public CellConvectiveOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    AdvDiffCenteredConvectiveOperator(std::string object_name,
                                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                      ConvectiveDifferencingType difference_form,
                                      std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> bc_coefs);

    /*!
     * \brief Default destructor.
     */
    ~AdvDiffCenteredConvectiveOperator() = default;

    /*!
     * \brief Static function to construct an AdvDiffCenteredConvectiveOperator.
     */
    static SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocate_operator(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      ConvectiveDifferencingType difference_form,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    {
        return new AdvDiffCenteredConvectiveOperator(object_name, Q_var, input_db, difference_form, bc_coefs);
    } // allocate_operator

    /*!
     * \brief Interpolate a cell-centered field Q to a face-centered field q on a single grid patch.
     */
    void interpolateToFaceOnPatch(SAMRAI::pdat::FaceData<NDIM, double>& q_interp_data,
                                  const SAMRAI::pdat::CellData<NDIM, double>& Q_cell_data,
                                  const SAMRAI::pdat::FaceData<NDIM, double>& u_data,
                                  const SAMRAI::hier::Patch<NDIM>& patch) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffCenteredConvectiveOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffCenteredConvectiveOperator(const AdvDiffCenteredConvectiveOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffCenteredConvectiveOperator& operator=(const AdvDiffCenteredConvectiveOperator& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_AdvDiffCenteredConvectiveOperator
