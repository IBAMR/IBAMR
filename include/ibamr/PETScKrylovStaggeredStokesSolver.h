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

#ifndef included_IBAMR_PETScKrylovStaggeredStokesSolver
#define included_IBAMR_PETScKrylovStaggeredStokesSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/KrylovLinearSolverStaggeredStokesSolverInterface.h"

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

namespace IBAMR
{
/*!
 * \brief Class PETScKrylovStaggeredStokesSolver is an extension of class
 * PETScKrylovLinearSolver that provides an implementation of the
 * StaggeredStokesSolver interface.
 */
class PETScKrylovStaggeredStokesSolver : public IBTK::PETScKrylovLinearSolver,
                                         public KrylovLinearSolverStaggeredStokesSolverInterface
{
public:
    /*!
     * \brief Constructor.
     */
    PETScKrylovStaggeredStokesSolver(const std::string& object_name,
                                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                     const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~PETScKrylovStaggeredStokesSolver() = default;

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScKrylovStaggeredStokesSolver() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScKrylovStaggeredStokesSolver(const PETScKrylovStaggeredStokesSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScKrylovStaggeredStokesSolver& operator=(const PETScKrylovStaggeredStokesSolver& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_PETScKrylovStaggeredStokesSolver
