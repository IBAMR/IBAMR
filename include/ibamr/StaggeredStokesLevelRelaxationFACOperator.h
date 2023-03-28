// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_StaggeredStokesLevelRelaxationFACOperator
#define included_IBAMR_StaggeredStokesLevelRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"

#include "ibtk/FACPreconditionerStrategy.h"

#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

#include <array>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoxList;
} // namespace hier
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
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
 * \brief Class StaggeredStokesLevelRelaxationFACOperator is a concrete
 * StaggeredStokesFACPreconditionerStrategy implementing a level relaxation
 * smoother for use as a multigrid preconditioner.
 */
class StaggeredStokesLevelRelaxationFACOperator : public StaggeredStokesFACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    StaggeredStokesLevelRelaxationFACOperator(const std::string& object_name,
                                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                              const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesLevelRelaxationFACOperator();

    /*!
     * \brief Static function to construct a StaggeredStokesFACPreconditioner with a
     * StaggeredStokesLevelRelaxationFACOperator FAC strategy.
     */
    static SAMRAI::tbox::Pointer<StaggeredStokesSolver>
    allocate_solver(const std::string& object_name,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    const std::string& default_options_prefix)
    {
        SAMRAI::tbox::Pointer<StaggeredStokesFACPreconditionerStrategy> fac_operator =
            new StaggeredStokesLevelRelaxationFACOperator(
                object_name + "::StaggeredStokesLevelRelaxationFACOperator", input_db, default_options_prefix);
        return new StaggeredStokesFACPreconditioner(object_name, fac_operator, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * \name Implementation of FACPreconditionerStrategy interface.
     */
    //\{

    /*!
     * \brief Specify the smoother type.
     */
    void setSmootherType(const std::string& level_solver_type);

    /*!
     * \brief Perform a given number of relaxations on the error.
     *
     * \param error error vector
     * \param residual residual vector
     * \param level_num level number
     * \param num_sweeps number of sweeps to perform
     * \param performing_pre_sweeps boolean value that is true when pre-smoothing sweeps are
     *being
     *performed
     * \param performing_post_sweeps boolean value that is true when post-smoothing sweeps are
     *being
     *performed
     */
    void smoothError(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                     const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                     int level_num,
                     int num_sweeps,
                     bool performing_pre_sweeps,
                     bool performing_post_sweeps) override;

    //\}

protected:
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    void initializeOperatorStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                            const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                                            int coarsest_reset_ln,
                                            int finest_reset_ln) override;

    /*!
     * \brief Remove implementation-specific hierarchy-dependent data.
     */
    void deallocateOperatorStateSpecialized(int coarsest_reset_ln, int finest_reset_ln) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesLevelRelaxationFACOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesLevelRelaxationFACOperator(const StaggeredStokesLevelRelaxationFACOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesLevelRelaxationFACOperator&
    operator=(const StaggeredStokesLevelRelaxationFACOperator& that) = delete;

    /*
     * Level solvers and solver parameters.
     */
    std::string d_level_solver_type = StaggeredStokesSolverManager::PETSC_LEVEL_SOLVER,
                d_level_solver_default_options_prefix;
    double d_level_solver_abs_residual_tol = 1.0e-50, d_level_solver_rel_residual_tol = 1.0e-5;
    int d_level_solver_max_iterations = 1;
    std::vector<SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesSolver> > d_level_solvers;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_level_solver_db;

    /*
     * Mappings from patch indices to patch operators.
     */
    std::vector<std::vector<std::array<SAMRAI::hier::BoxList<NDIM>, NDIM> > > d_patch_side_bc_box_overlap;
    std::vector<std::vector<SAMRAI::hier::BoxList<NDIM> > > d_patch_cell_bc_box_overlap;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_StaggeredStokesLevelRelaxationFACOperator
