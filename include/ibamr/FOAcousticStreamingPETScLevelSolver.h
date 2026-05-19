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

#ifndef included_IBAMR_FOAcousticStreamingPETScLevelSolver
#define included_IBAMR_FOAcousticStreamingPETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/PETScLevelSolver.h"
#include "ibtk/ibtk_utilities.h"

#include "IntVector.h"
#include "RefineSchedule.h"
#include "SideVariable.h"
#include "VariableContext.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "petscvec.h"

#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class FOAcousticStreamingPETScLevelSolver is a parent PETScLevelSolver for
 * solving the first order acoustic streaming equations on a \em single SAMRAI::hier::PatchLevel
 * using <A HREF="http://www.mcs.anl.gov/petsc/petsc-as">PETSc</A>.
 *
 * This solver class uses the PETSc library to solve linear equations for
 * both cell- and side-centered velocity and cell-centered pressure. Both velocity and pressure
 * fields are complex valued functions, and have two components (real and imaginary) each.
 * The discretization is second-order accurate.
 */
class FOAcousticStreamingPETScLevelSolver : public IBTK::PETScLevelSolver
{
public:
    /*!
     * \brief Destructor.
     */
    ~FOAcousticStreamingPETScLevelSolver();

    /*!
     * \brief Set mass density patch data index.
     */
    virtual void setMassDensityPatchDataIndex(int rho_idx)
    {
        d_rho_idx = rho_idx;
        return;
    } // registerMassDensityPatchDataIndex

    /*!
     * \brief Set shear/first coefficient of viscosity patch data index.
     */
    virtual void setShearViscosityPatchDataIndex(int mu_idx)
    {
        d_mu_idx = mu_idx;
        return;
    } // setShearViscosityPatchDataIndex

    /*!
     * \brief Set bulk/second coefficient of viscosity patch data index.
     */
    virtual void setBulkViscosityPatchDataIndex(int lambda_idx)
    {
        d_lambda_idx = lambda_idx;
        return;
    } // setBulkViscosityPatchDataIndex

    virtual void setViscosityInterpolationType(const IBTK::VCInterpType mu_interp_type)
    {
        d_mu_interp_type = mu_interp_type;
        return;
    } // setViscosityInterpolationType

    /*!
     * \brief Set Brinkman penalization coefficient of velocity patch data index.
     */
    virtual void setBrinkmanPenalizationPatchDataIndex(int chi_idx)
    {
        d_chi_idx = chi_idx;
        return;
    } // setBrinkmanPenalizationPatchDataIndex

    /*!
     * \brief Set acoustic angular frequency
     */
    virtual void setAcousticAngularFrequency(double omega)
    {
        d_omega = omega;
        return;
    } // setAcousticAngularFrequency

    /*!
     * \brief Set the speed of sound in the fluid medium
     */
    virtual void setSoundSpeed(double c0)
    {
        d_sound_speed = c0;
        return;
    } // setSoundSpeed

    /*!
     * \brief Set boundary conditions for the velocity components
     */
    virtual void setBoundaryConditionCoefficients(
        const std::array<std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>, 2>& U_bc_coefs)
    {
        d_U_bc_coefs = U_bc_coefs;
        d_P_bc_coefs = { nullptr, nullptr };
        return;
    } // setBoundaryConditionCoefficients

    /*!
     * \brief Set boundary conditions for the velocity and pressure components
     */
    virtual void setBoundaryConditionCoefficients(
        const std::array<std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>, 2>& U_bc_coefs,
        const std::array<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*, 2>& P_bc_coefs)
    {
        d_U_bc_coefs = U_bc_coefs;
        d_P_bc_coefs = P_bc_coefs;
        return;
    } // setBoundaryConditionCoefficients

protected:
    /*!
     * \brief Constructor.
     */
    FOAcousticStreamingPETScLevelSolver(const std::string& object_name,
                                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > u_dof_idx_var,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > p_dof_idx_var,
                                        const std::string& default_options_prefix);

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     */
    void initializeSolverStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                          const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverStateSpecialized().
     */
    void deallocateSolverStateSpecialized() override;

    /*!
     * \brief Copy a generic vector to the PETSc representation.
     */
    void copyToPETScVec(Vec& petsc_x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x) override;

    /*!
     * \brief Copy a generic vector from the PETSc representation.
     */
    void copyFromPETScVec(Vec& petsc_x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x) override;

    /*!
     * \name PETSc objects.
     */
    //\{
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    std::vector<int> d_num_dofs_per_proc;
    int d_u_dof_index_idx = IBTK::invalid_index, d_p_dof_index_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_u_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, int> > d_p_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > d_data_synch_sched, d_ghost_fill_sched;
    //\}

    /*!
     * \name Parameters specific to the first order acoustic streaming system.
     */
    double d_omega = std::numeric_limits<double>::signaling_NaN(),
           d_sound_speed = std::numeric_limits<double>::signaling_NaN();

    /*!
     * \name Patch data indices for material parameters.
     */
    int d_rho_idx = IBTK::invalid_index, d_lambda_idx = IBTK::invalid_index, d_mu_idx = IBTK::invalid_index;

    /*!
     * \name Patch data index for the Brinkman penalization coefficient.
     */
    int d_chi_idx = IBTK::invalid_index;

    /*!
     * Velocity boundary conditions for the real and imaginary components
     */
    std::array<std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>, 2> d_U_bc_coefs;

    /*!
     * Pressure boundary conditions for the real and imaginary components
     */
    std::array<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*, 2> d_P_bc_coefs;

    /*!
     * Interpolation type for shear viscosity
     */
    IBTK::VCInterpType d_mu_interp_type = IBTK::VC_HARMONIC_INTERP;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    FOAcousticStreamingPETScLevelSolver() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FOAcousticStreamingPETScLevelSolver(const FOAcousticStreamingPETScLevelSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FOAcousticStreamingPETScLevelSolver& operator=(const FOAcousticStreamingPETScLevelSolver& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_FOAcousticStreamingPETScLevelSolver
