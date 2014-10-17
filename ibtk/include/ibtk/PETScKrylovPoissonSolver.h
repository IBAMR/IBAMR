// Filename: PETScKrylovPoissonSolver.h
// Created on 20 Nov 2012 by Boyce Griffith
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

#ifndef included_PETScKrylovPoissonSolver
#define included_PETScKrylovPoissonSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "ibtk/KrylovLinearSolverPoissonSolverInterface.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "tbox/Pointer.h"

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
    PETScKrylovPoissonSolver(const std::string& object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~PETScKrylovPoissonSolver();

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScKrylovPoissonSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScKrylovPoissonSolver(const PETScKrylovPoissonSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScKrylovPoissonSolver& operator=(const PETScKrylovPoissonSolver& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScKrylovPoissonSolver
