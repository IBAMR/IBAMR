// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_NewtonKrylovSolverManager
#define included_IBTK_NewtonKrylovSolverManager

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/NewtonKrylovSolver.h"

#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class NewtonKrylovSolverManager is a singleton manager class to
 * provide access to generic NewtonKrylovSolver implementations.
 */
class NewtonKrylovSolverManager
{
public:
    /*!
     * Key to use for "undefined" solver types.
     */
    static const std::string UNDEFINED;

    /*!
     * Default Newton-Krylov solver types automatically provided by the manager
     * class.
     */
    static const std::string DEFAULT;
    static const std::string PETSC;

    /*!
     * Return a pointer to the instance of the solver manager.  Access to
     * NewtonKrylovSolverManager objects is mediated by the getManager()
     * function.
     *
     * \return A pointer to the solver manager instance.
     */
    static NewtonKrylovSolverManager* getManager();

    /*!
     * Deallocate the NewtonKrylovSolverManager instance.
     *
     * It is not necessary to call this function at program termination since it
     * is automatically called by the ShutdownRegistry class.
     */
    static void freeManager();

    /*!
     * Allocate a new NewtonKrylovSolver object of the specified type.
     */
    SAMRAI::tbox::Pointer<NewtonKrylovSolver>
    allocateSolver(const std::string& solver_type,
                   const std::string& solver_object_name,
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> solver_input_db,
                   const std::string& solver_default_options_prefix) const;

    /*!
     * Typedef for functions to construct NewtonKrylovSolvers.
     */
    using SolverMaker =
        SAMRAI::tbox::Pointer<NewtonKrylovSolver> (*)(const std::string& solver_object_name,
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
    NewtonKrylovSolverManager();

    /*!
     * \brief Destructor.
     */
    ~NewtonKrylovSolverManager() = default;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    NewtonKrylovSolverManager(const NewtonKrylovSolverManager& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    NewtonKrylovSolverManager& operator=(const NewtonKrylovSolverManager& that) = delete;

    /*!
     * Static data members used to control access to and destruction of
     * singleton data manager instance.
     */
    static NewtonKrylovSolverManager* s_solver_manager_instance;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    /*!
     * Mapping from solver type names to solver maker functions.
     */
    std::map<std::string, SolverMaker> d_solver_maker_map;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_NewtonKrylovSolverManager
