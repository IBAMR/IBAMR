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

#ifndef included_EnthalpyCorrectionBcCoefs
#define included_EnthalpyCorrectionBcCoefs

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Patch.h>
#include <RobinBcCoefStrategy.h>
#include <Variable.h>

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Robin boundary coefficients for an enthalpy Newton correction.
 *
 * For a temperature boundary condition a_T T + b_T dT/dn = g_T, this class
 * supplies a_h = a_T beta, b_h = b_T beta, and
 * g_h = g_T - a_T T^m - b_T dT^m/dn, where beta = dT/dh.
 *
 * The composition and constitutive derivative are lagged during each linear
 * solve. The term b_T (d beta/dn) delta_h is neglected.
 *
 * Temperature coefficients are evaluated by the parent muParser object or
 * an optional externally supplied temperature Robin boundary object.
 *
 * \note The cell-centered temperature and derivative data must be allocated,
 * ghost-filled, and frozen throughout each correction solve. In particular,
 * the temperature index must not be the scratch index used for delta_h.
 */
class EnthalpyCorrectionBcCoefs : public SAMRAI::solv::RobinBcCoefStrategy<NDIM>

{
    using Array = SAMRAI::pdat::ArrayData<NDIM, double>;
    using ArrayPointer = SAMRAI::tbox::Pointer<Array>;
    using VariablePointer = SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>;

public:
    /*!
     * \brief Construct using the temperature Robin coefficient database.
     */
    EnthalpyCorrectionBcCoefs(const std::string& object_name);

    /*!
     * \brief Destructor.
     */
    ~EnthalpyCorrectionBcCoefs() override = default;

    /*!
     * \brief Set the cell-centered temperature and dT/dh patch-data indices.
     */
    void setTemperatureAndDerivativePatchDataIndices(int T_idx, int dT_dh_idx);

    /*!
     * \brief Delegate temperature coefficients to an existing boundary condition object.
     *
     * The pointer is non-owning. The supplied object must outlive all boundary evaluations.
     */
    void setTemperatureBcCoef(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* temperature_bc_coef);

    /*!
     * \brief Fill the Robin coefficients for delta_h.
     *
     * Any output array may be null. Boundary temperature and beta are sampled
     * by averaging the adjacent interior and ghost values. The outward normal
     * temperature derivative is (T_ghost - T_interior)/dx.
     */
    void setBcCoefs(ArrayPointer& acoef_data,
                    ArrayPointer& bcoef_data,
                    ArrayPointer& gcoef_data,
                    const VariablePointer& variable,
                    const SAMRAI::hier::Patch<NDIM>& patch,
                    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                    double fill_time = 0.0) const override;

    /*!
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     */
    SAMRAI::hier::IntVector<NDIM> numberOfExtensionsFillable() const override;

private:
    EnthalpyCorrectionBcCoefs() = delete;
    EnthalpyCorrectionBcCoefs(const EnthalpyCorrectionBcCoefs&) = delete;
    EnthalpyCorrectionBcCoefs& operator=(const EnthalpyCorrectionBcCoefs&) = delete;

    void fillConvertedCoefficient(ArrayPointer& acoef_data,
                                  ArrayPointer& bcoef_data,
                                  ArrayPointer& gcoef_data,
                                  const VariablePointer& variable,
                                  const SAMRAI::hier::Patch<NDIM>& patch,
                                  const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                                  double fill_time) const;

    std::string d_object_name;
    int d_T_idx = IBTK::invalid_index;
    int d_dT_dh_idx = IBTK::invalid_index; // dT_dh patch-data index, cell-centered

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_temperature_bc_coef = nullptr;
};
} // namespace IBAMR

#endif // included_EnthalpyCorrectionBcCoefs
