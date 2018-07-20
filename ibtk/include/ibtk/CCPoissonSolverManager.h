// Filename: CCPoissonSolverManager.h
// Created on 13 Aug 2012 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#ifndef included_IBTK_CCPoissonSolverManager
#define included_IBTK_CCPoissonSolverManager

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <string>

#include "ibtk/PoissonSolver.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CCPoissonSolverManager is a singleton manager class to provide
 * access to generic cell-centered PoissonSolver implementations.
 */
class CCPoissonSolverManager
{
public:
    /*!
     * Key to use for "undefined" solver types.
     */
    static const std::string UNDEFINED;

    /*!
     * Default Krylov solver types automatically provided by the manager class.
     */
    static const std::string DEFAULT_KRYLOV_SOLVER;
    static const std::string PETSC_KRYLOV_SOLVER;

    /*!
     * Default FAC preconditioner types automatically provided by the manager
     * class.
     */
    static const std::string DEFAULT_FAC_PRECONDITIONER;
    static const std::string BOX_RELAXATION_FAC_PRECONDITIONER;
    static const std::string LEVEL_RELAXATION_FAC_PRECONDITIONER;
    static const std::string POINT_RELAXATION_FAC_PRECONDITIONER;

    /*!
     * Default level solver types automatically provided by the manager class.
     */
    static const std::string DEFAULT_LEVEL_SOLVER;
    static const std::string HYPRE_LEVEL_SOLVER;
    static const std::string PETSC_LEVEL_SOLVER;

    /*!
     * Return a pointer to the instance of the solver manager.  Access to
     * CCPoissonSolverManager objects is mediated by the getManager()
     * function.
     *
     * \return A pointer to the solver manager instance.
     */
    static CCPoissonSolverManager* getManager();

    /*!
     * Deallocate the CCPoissonSolverManager instance.
     *
     * It is not necessary to call this function at program termination since it
     * is automatically called by the ShutdownRegistry class.
     */
    static void freeManager();

    /*!
     * Allocate a new CCPoissonSolver object of the specified type.
     */
    SAMRAI::tbox::Pointer<PoissonSolver> allocateSolver(const std::string& solver_type,
                                                        const std::string& solver_object_name,
                                                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                                                        const std::string& solver_default_options_prefix) const;

    /*!
     * Allocate a new CCPoissonSolver object of the specified type with a
     * preconditioner of a specified type.
     *
     * \note The preconditioner settings are used only when the parent solver
     * is a KrylovLinearSolver.
     */
    SAMRAI::tbox::Pointer<PoissonSolver>
    allocateSolver(const std::string& solver_type,
                   const std::string& solver_object_name,
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                   const std::string& solver_default_options_prefix,
                   const std::string& precond_type,
                   const std::string& precond_object_name,
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> precond_input_db,
                   const std::string& precond_default_options_prefix,
                   const std::string& sub_precond_type = "",
                   const std::string& sub_precond_object_name = "",
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> sub_precond_input_db =
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(),
                   const std::string& sub_precond_default_options_prefix = "") const;

    /*!
     * Typedef for functions to construct cell-centered PoissonSolvers.
     */
    typedef SAMRAI::tbox::Pointer<PoissonSolver> (*SolverMaker)(
        const std::string& solver_object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
        const std::string& solver_default_options_prefix);

    /*!
     * Register a solver factory function with the solver manager class.
     */
    void registerSolverFactoryFunction(const std::string& solver_type, SolverMaker solver_maker);

protected:
    /*!
     * \brief Default constructor.
     */
    CCPoissonSolverManager();

    /*!
     * \brief Destructor.
     */
    ~CCPoissonSolverManager();

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CCPoissonSolverManager(const CCPoissonSolverManager& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CCPoissonSolverManager& operator=(const CCPoissonSolverManager& that) = delete;

    /*!
     * Static data members used to control access to and destruction of
     * singleton data manager instance.
     */
    static CCPoissonSolverManager* s_solver_manager_instance;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    /*!
     * Mapping from solver type names to solver maker functions.
     */
    std::map<std::string, SolverMaker> d_solver_maker_map;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CCPoissonSolverManager
