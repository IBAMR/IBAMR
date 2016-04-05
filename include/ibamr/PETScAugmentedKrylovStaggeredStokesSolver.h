// Filename: PETScAugmentedKrylovStaggeredStokesSolver.h
// Created on 04 Apr 2016 by Boyce Griffith
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

#ifndef included_PETScAugmentedKrylovStaggeredStokesSolver
#define included_PETScAugmentedKrylovStaggeredStokesSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "ibamr/KrylovLinearSolverStaggeredStokesSolverInterface.h"
#include "ibtk/PETScAugmentedKrylovLinearSolver.h"
#include "tbox/Pointer.h"

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
 * \brief Class PETScAugmentedKrylovStaggeredStokesSolver is an extension of class
 * PETScAugmentedKrylovLinearSolver that provides an implementation of the
 * StaggeredStokesSolver interface.
 */
class PETScAugmentedKrylovStaggeredStokesSolver : public IBTK::PETScAugmentedKrylovLinearSolver,
                                                  public KrylovLinearSolverStaggeredStokesSolverInterface
{
public:
    /*!
     * \brief Constructor.
     */
    PETScAugmentedKrylovStaggeredStokesSolver(const std::string& object_name,
                                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                              const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~PETScAugmentedKrylovStaggeredStokesSolver();

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScAugmentedKrylovStaggeredStokesSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScAugmentedKrylovStaggeredStokesSolver(const PETScAugmentedKrylovStaggeredStokesSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScAugmentedKrylovStaggeredStokesSolver& operator=(const PETScAugmentedKrylovStaggeredStokesSolver& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScAugmentedKrylovStaggeredStokesSolver
