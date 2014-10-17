// Filename: StokesBcCoefStrategy.h
// Created on 04 Sep 2012 by Boyce Griffith
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

#ifndef included_StokesBcCoefStrategy
#define included_StokesBcCoefStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/ibamr_enums.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"

namespace IBAMR
{
class StokesSpecifications;
} // namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StokesBcCoefStrategy is a subclass of the abstract base class
 * IBTK::ExtendedRobinBcCoefStrategy to allow for specialization needed to treat
 * physical boundary conditions for coupled Stokes solvers.
 */
class StokesBcCoefStrategy : public IBTK::ExtendedRobinBcCoefStrategy
{
public:
    /*!
     * \brief Empty default constructor.
     */
    StokesBcCoefStrategy();

    /*!
     * \brief Empty destructor.
     */
    ~StokesBcCoefStrategy();

    /*!
     * \brief Set the StokesSpecifications object used by this boundary condition
     * specification object.
     *
     * \param problem_coefs   Problem coefficients
     */
    virtual void setStokesSpecifications(const StokesSpecifications* problem_coefs);

    /*!
     * \brief Set the target velocity data index to use when setting physical
     * boundary conditions and the time at which it is defined.
     */
    virtual void setTargetVelocityPatchDataIndex(int u_target_data_idx);

    /*!
     * \brief Clear the target velocity data index used when setting physical
     * boundary conditions.
     */
    virtual void clearTargetVelocityPatchDataIndex();

    /*!
     * \brief Set the target pressure data index to use when setting physical
     * boundary conditions and the time at which it is defined.
     */
    virtual void setTargetPressurePatchDataIndex(int u_target_data_idx);

    /*!
     * \brief Clear the target pressure data index used when setting physical
     * boundary conditions.
     */
    virtual void clearTargetPressurePatchDataIndex();

    /*!
     * \brief Set the type of traction boundary conditions.  Supported options
     * are: TRACTION_BOUNDARY_CONDITIONS and
     * PSEUDO_TRACTION_BOUNDARY_CONDITIONS.
     *
     * The default is TRACTION_BOUNDARY_CONDITIONS.
     */
    virtual void setTractionBcType(TractionBcType bc_type);

    /*!
     * \brief Get the type of traction boundary conditions.
     */
    virtual TractionBcType getTractionBcType() const;

protected:
    /*
     * Problem coefficients.
     */
    const StokesSpecifications* d_problem_coefs;

    /*!
     * Patch data indices.
     */
    int d_u_target_data_idx, d_p_target_data_idx;

    /*
     * The type of traction boundary conditions.  Supported options are:
     * TRACTION_BOUNDARY_CONDITIONS and PSEUDO_TRACTION_BOUNDARY_CONDITIONS.
     */
    TractionBcType d_traction_bc_type;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StokesBcCoefStrategy(const StokesBcCoefStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StokesBcCoefStrategy& operator=(const StokesBcCoefStrategy& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StokesBcCoefStrategy
