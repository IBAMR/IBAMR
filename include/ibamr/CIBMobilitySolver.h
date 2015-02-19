// Filename: CIBMobilitySolver.h
// Created on 19 Feb 2015 by Amneet Bhalla and Bakytzhan Kallemov
//
// Copyright (c) 2002-2015, Amneet Bhalla, Bakytzhan Kallemov,
// and Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_CIBIBMobilitySolver
#define included_CIBIBMobilitySolver

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <vector>

#include "petscksp.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "tbox/Database.h"

namespace IBAMR
{
class INSStaggeredHierarchyIntegrator;
class CIBStrategy;
class DenseMobilitySolver;
class KrylovMobilitySolver;
}// namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBAMR
{

/*!
 * \brief Class CIBMobilitySolver solves for the constraint forces 
 * \f$ \vec{\lambda}\f$ required to impose the rigid velocity as defined by 
 * the mobility sub-problem. Specifically, the class solves the matrix-equation
 * 
 *  \f{eqnarray*}{
 *      J L^{-1} S \vec{\lambda} = M \vec{\lambda} = \vec{U}.
 * \f}
 *
 * Here, \f$ J \f$ is the interpolation operator, \f$ S \f$ is the spreading 
 * operator, \f$ L^{-1} \f$ is inverse of incompressible Stokes operator, and
 * \f$ \vec{U} \f$ is the desired velocity of the structure(s).
 *
 * This class employs direct solver for the approximate mobility sub-problem.
 * The approximate mobility matrix is intended to be used in the preconditioning
 * step of the overall constraint solver as defined in IBAMR::CIBSaddlePointSolver 
 * class. The class also supports (an unpreconditioned) Krylov mobility solver for
 * small test problems.
 */

enum MOBILITY_SOLVER_TYPE
{
DIRECT = 0,
KRYLOV,
SOL_TYPE_UNKNOWN = -1
};

enum MOBILITY_SOLVER_SUBTYPE
{
DENSE = 0,
BLOCK_DIAGONAL,
KRYLOV_REGULAR,
KRYLOV_FMM,
SOL_SUBTYPE_UNKNOWN = -1
};
	
class CIBMobilitySolver
    : public SAMRAI::tbox::DescribedClass
{

/////////////////////////////// PUBLIC ///////////////////////////////////////

public:
	
    /*!
     * \brief The only constructor of this class.
     */
    CIBMobilitySolver(
		const std::string& object_name,
		SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
		SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
		SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> cib_strategy);

    /*!
     * \brief Destructor for this class.
     */
    ~CIBMobilitySolver();

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void
    setSolutionTime(
		const double solution_time);

    /*!
     * \brief Set the time interval of integration.
     */
    void
    setTimeInterval(
		const double current_time,
		const double new_time);
   
    /*!
     * \brief Solves the mobility problem.
	 *
	 * \param x Vec storing the Lagrange multiplier
	 *
	 * \param b Vec storing the desired velocity
	 *
	 * \return \p true if the solver converged to the specified tolerances, \p
	 * false otherwise
     */
    bool
    solveSystem(
        Vec x, 
		Vec b);

    /*!
     * \brief Initialize the mobility solver.
	 *
	 * \param x Vec storing the Lagrange multiplier
	 *
	 * \param b Vec storing the desired velocity
     *   
     */
    void
    initializeSolverState(
		Vec x,
		Vec b);

    /*!
     * \brief Deallocate the mobility solver.
     *   
     */
    void
    deallocateSolverState();
	
/////////////////////////////// PRIVATE //////////////////////////////////////

private:
	
	/*!
	 * \brief Get various options from input db.
	 */
	void
	getFromInput(
		SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
	
    /*!
     * \brief Initialize the direct solver(s).
     */
    void
    initializeDirectSolver();
    
    /*!
     * \brief Collect the solution from various solvers to return it to 
	 * CIBSaddlePoint solver.
     */
    void
    returnSolution(
        Vec x); 

    //Name of this object.
    std::string d_object_name;
    
    // Pointers.
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_input_db;
    SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator> d_ins_integrator;
    SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> d_cib_strategy;
	std::vector<SAMRAI::tbox::Pointer<IBAMR::DenseMobilitySolver> > d_dense_solvers;
	SAMRAI::tbox::Pointer<IBAMR::KrylovMobilitySolver> d_krylov_solver;
    std::vector<Vec> d_nodes_x, d_nodes_w;

    // Other parameters.
    int d_num_rigid_parts;
    double d_solution_time, d_current_time, d_new_time, d_dt;
    bool d_recompute_mobility_matrix;         // recompute mobility each iteration
    bool d_is_initialized;                    // flag for initialization of mobility matrix
	
	// Type of mobility solver to be used
    MOBILITY_SOLVER_TYPE    d_mobility_solver_type;
	MOBILITY_SOLVER_SUBTYPE d_mobility_solver_subtype;
	
};// CIBMobilitySolver

}// IBAMR

#endif // #ifndef included_CIBMobilitySolver
