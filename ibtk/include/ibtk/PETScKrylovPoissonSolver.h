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

#ifndef included_IBTK_PETScKrylovPoissonSolver
#define included_IBTK_PETScKrylovPoissonSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/KrylovLinearSolverPoissonSolverInterface.h"
#include "ibtk/PETScKrylovLinearSolver.h"

#include "tbox/Pointer.h"

#include <string>

namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScKrylovPoissonSolver is an extension of class
 * PETScKrylovLinearSolver that provides an implementation of the PoissonSolver
 * interface.
 */
class PETScKrylovPoissonSolver : public PETScKrylovLinearSolver, public KrylovLinearSolverPoissonSolverInterface
{
public:
    /*!
     * \brief Constructor.
     */
    PETScKrylovPoissonSolver(std::string object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             std::string default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~PETScKrylovPoissonSolver() = default;

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScKrylovPoissonSolver() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScKrylovPoissonSolver(const PETScKrylovPoissonSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScKrylovPoissonSolver& operator=(const PETScKrylovPoissonSolver& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PETScKrylovPoissonSolver
