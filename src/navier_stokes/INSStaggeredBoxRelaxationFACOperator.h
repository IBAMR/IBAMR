// Filename: INSStaggeredBoxRelaxationFACOperator.h
// Created on 11 Jun 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
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

#ifndef included_INSStaggeredBoxRelaxationFACOperator
#define included_INSStaggeredBoxRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscksp.h>

// IBAMR INCLUDES
#include <ibamr/INSStaggeredFACPreconditionerStrategy.h>
#include <ibamr/INSProblemCoefs.h>

// SAMRAI INCLUDES
#include <LocationIndexRobinBcCoefs.h>

// C++ STDLIB INCLUDES
#include <map>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredBoxRelaxationFACOperator is a concrete
 * INSStaggeredFACPreconditionerStrategy implementing a box relaxation
 * (Vanka-type) smoother for use as a multigrid preconditioner.
*/
class INSStaggeredBoxRelaxationFACOperator
    : public INSStaggeredFACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    INSStaggeredBoxRelaxationFACOperator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db=NULL);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredBoxRelaxationFACOperator();

    /*!
     * \name Functions for specifying the problem coefficients.
     */
    //\{

    /*!
     * \brief Set the INSProblemCoefs object and timestep size used to specify
     * the coefficients for the time-dependent incompressible Stokes operator.
     */
    void
    setProblemCoefficients(
        const INSProblemCoefs& problem_coefs);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a U_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param U_bc_coefs  Vector of pointers to objects that can set the Robin boundary condition coefficients for the velocity
     * \param P_bc_coef   Pointer to object that can set the Robin boundary condition coefficients for the pressure
     */
    void
    setPhysicalBcCoefs(
        const blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM>& U_bc_coefs,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef);

    //\}

    /*!
     * \name Implementation of FACPreconditionerStrategy interface.
     */
    //\{

    /*!
     * \brief Perform a given number of relaxations on the error.
     *
     * \param error error vector
     * \param residual residual vector
     * \param level_num level number
     * \param num_sweeps number of sweeps to perform
     * \param performing_pre_sweeps boolean value that is true when pre-smoothing sweeps are being performed
     * \param performing_post_sweeps boolean value that is true when post-smoothing sweeps are being performed
     */
    void
    smoothError(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        int level_num,
        int num_sweeps,
        bool performing_pre_sweeps,
        bool performing_post_sweeps);

    /*!
     * \brief Solve the residual equation Ae=r on the coarsest level of the
     * patch hierarchy.
     *
     * \param error error vector
     * \param residual residual vector
     * \param coarsest_ln coarsest level number
     */
    bool
    solveCoarsestLevel(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        int coarsest_ln);

    /*!
     * \brief Compute the composite-grid residual on the specified range of
     * levels of the patch hierarchy.
     */
    void
    computeResidual(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs,
        int coarsest_level_num,
        int finest_level_num);

    //\}

protected:
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    void
    initializeOperatorStateSpecialized(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs,
        int coarsest_reset_ln,
        int finest_reset_ln);

    /*!
     * \brief Remove implementation-specific hierarchy-dependent data.
     */
    void
    deallocateOperatorStateSpecialized(
        int coarsest_reset_ln,
        int finest_reset_ln);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredBoxRelaxationFACOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredBoxRelaxationFACOperator(
        const INSStaggeredBoxRelaxationFACOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredBoxRelaxationFACOperator& operator=(
        const INSStaggeredBoxRelaxationFACOperator& that);

    /*
     * Problem coefficient specifications.
     */
    INSProblemCoefs d_problem_coefs;

    /*
     * \name Boundary condition handling objects.
     */
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* const d_default_U_bc_coef;
    blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM> d_U_bc_coefs;

    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* const d_default_P_bc_coef;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef;

    /*
     * Box operator data.
     */
    std::vector<Mat> d_box_op;
    std::vector<Vec> d_box_e, d_box_r;
    std::vector<KSP> d_box_ksp;

    /*
     * Mappings from patch indices to patch operators.
     */
    std::vector<std::vector<blitz::TinyVector<SAMRAI::hier::BoxList<NDIM>,NDIM> > > d_patch_side_bc_box_overlap;
    std::vector<std::vector<SAMRAI::hier::BoxList<NDIM> > > d_patch_cell_bc_box_overlap;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/INSStaggeredBoxRelaxationFACOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredBoxRelaxationFACOperator
