// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
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

#ifndef included_IBTK_ProblemSpecification
#define included_IBTK_ProblemSpecification

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBTK
{

/*!
 * @brief An empty structure holding problem specifications for operators and solvers.
 *
 */
struct ProblemSpecification
{
};

/*!
 * \brief Class VCViscousOpSpec is a light-weight structure holding
 * values of VCViscousOp coefficients.
 */
struct VCViscousOpSpec : public ProblemSpecification
{
    VCViscousOpSpec()
    {
        return;
    } // VCViscousOpSpec

    /*
     * Reset the state of the problem specification object
     */
    void reset()
    {
        d_C_idx = IBTK::invalid_index;
        d_D_idx = IBTK::invalid_index;
        d_D_const = std::numeric_limits<double>::signaling_NaN();
        d_C_const = std::numeric_limits<double>::signaling_NaN();
        d_D_is_const = false;
        d_C_is_const = false;
    } // reset

    // C_idx is the intertial and/or Brinkman penalization coefficient
    // D_idx is the shear viscosity coefficient
    int d_C_idx = IBTK::invalid_index, d_D_idx = IBTK::invalid_index;

    double d_D_const = std::numeric_limits<double>::signaling_NaN();
    double d_C_const = std::numeric_limits<double>::signaling_NaN();

    bool d_D_is_const = false, d_C_is_const = false;

}; // VCViscousOpSpec

/*!
 * \brief Class VCViscousDilatationalOpSpec is a light-weight structure holding
 * values of VCViscousDilatationalOp coefficients.
 */
struct VCViscousDilatationalOpSpec : public VCViscousOpSpec
{
    VCViscousDilatationalOpSpec() : VCViscousOpSpec()
    {
        return;
    } // VCViscousDilatationalOpSpec

    /*
     * Reset the state of the problem specification object
     */
    void reset()
    {
        VCViscousOpSpec::reset();
        d_L_idx = IBTK::invalid_index;
        d_L_const = std::numeric_limits<double>::signaling_NaN();
        d_L_is_const = false;
    } // reset

    operator VCViscousOpSpec()
    {
        VCViscousOpSpec obj;
        obj.d_C_idx = this->d_C_idx;
        obj.d_D_idx = this->d_D_idx;
        obj.d_D_const = this->d_D_const;
        obj.d_C_const = this->d_C_const;
        obj.d_D_is_const = this->d_D_is_const;
        obj.d_C_is_const = this->d_C_is_const;

        return obj;
    } // operator VCViscousOpSpec

    // L_idx is the dilatational viscosity coefficient
    int d_L_idx = IBTK::invalid_index;
    double d_L_const = std::numeric_limits<double>::signaling_NaN();
    bool d_L_is_const = false;

}; // VCViscousDilatationalOpSpec

} // namespace IBTK
//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_ProblemSpecification
