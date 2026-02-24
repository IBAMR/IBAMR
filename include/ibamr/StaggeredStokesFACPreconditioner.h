// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#ifndef included_IBAMR_StaggeredStokesFACPreconditioner
#define included_IBAMR_StaggeredStokesFACPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StaggeredStokesSolver.h"

#include "ibtk/FACPreconditioner.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIDatabase.h"
#include "SAMRAIPointer.h"
#include "SAMRAIPoissonSpecifications.h"
#include "SAMRAIRobinBcCoefStrategy.h"

#include <string>
#include <vector>

namespace IBAMR
{
class StaggeredStokesPhysicalBoundaryHelper;
} // namespace IBAMR
namespace IBTK
{
class FACPreconditionerStrategy;
} // namespace IBTK
namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesFACPreconditioner is a IBTK::FACPreconditioner that has
 * been specialized for Stokes problems.
 */
class StaggeredStokesFACPreconditioner : public IBTK::FACPreconditioner, public StaggeredStokesSolver
{
public:
    /*!
     * Constructor.
     */
    StaggeredStokesFACPreconditioner(const std::string& object_name,
                                     SAMRAIPointer<IBTK::FACPreconditionerStrategy> fac_strategy,
                                     SAMRAIPointer<SAMRAIDatabase> input_db,
                                     const std::string& default_options_prefix);

    /*!
     * Destructor.
     */
    ~StaggeredStokesFACPreconditioner() = default;

    /*!
     * \brief Set the SAMRAIPoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    void setVelocityPoissonSpecifications(const SAMRAIPoissonSpecifications& U_problem_coefs) override;

    /*!
     * \brief Set if velocity and pressure have nullspace.
     */
    void setComponentsHaveNullSpace(const bool has_velocity_nullspace, const bool has_pressure_nullspace) override;

    /*!
     * \brief Set the SAMRAIRobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a U_bc_coefs may be nullptr.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be nullptr; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param U_bc_coefs  IBTK::Vector of pointers to objects that can set the Robin boundary
     *condition coefficients for the velocity
     * \param P_bc_coef   SAMRAIPointer to object that can set the Robin boundary condition
     *coefficients
     *for the pressure
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAIRobinBcCoefStrategy*>& U_bc_coefs,
                            SAMRAIRobinBcCoefStrategy* P_bc_coef) override;

    /*!
     * \brief Set the StokesSpecifications object and timestep size used to specify
     * the coefficients for the time-dependent incompressible Stokes operator.
     */
    void setPhysicalBoundaryHelper(SAMRAIPointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesFACPreconditioner() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesFACPreconditioner(const StaggeredStokesFACPreconditioner& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesFACPreconditioner& operator=(const StaggeredStokesFACPreconditioner& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesFACPreconditioner
