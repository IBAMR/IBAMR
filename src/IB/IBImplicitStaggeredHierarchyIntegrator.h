// Filename: IBImplicitStaggeredHierarchyIntegrator.h
// Created on 07 Apr 2012 by Boyce Griffith
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

#ifndef included_IBImplicitStaggeredHierarchyIntegrator
#define included_IBImplicitStaggeredHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSC INCLUDES
#include <petscsys.h>

// IBAMR INCLUDES
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

// IBTK INCLUDES
#include <ibtk/NewtonKrylovSolver.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBImplicitStaggeredHierarchyIntegrator is an implementation of a
 * formally second-order accurate, nonlinearly-implicit version of the immersed
 * boundary method.
 */
class IBImplicitStaggeredHierarchyIntegrator
    : public IBHierarchyIntegrator
{
public:
    /*!
     * The constructor for class IBImplicitStaggeredHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    IBImplicitStaggeredHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<IBStrategy> ib_method_ops,
        SAMRAI::tbox::Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
        bool register_for_restart=true);

    /*!
     * The destructor for class IBImplicitStaggeredHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
     */
    ~IBImplicitStaggeredHierarchyIntegrator();

    /*!
     * Returns the number of cycles to perform for the present time step.
     */
    int
    getNumberOfCycles();

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void
    preprocessIntegrateHierarchy(
        double current_time,
        double new_time,
        int num_cycles=1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void
    integrateHierarchy(
        double current_time,
        double new_time,
        int cycle_num=0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void
    postprocessIntegrateHierarchy(
        double current_time,
        double new_time,
        bool skip_synchronize_new_state_data,
        int num_cycles=1);

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void
    initializeHierarchyIntegrator(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

protected:
    /*!
     * Write out specialized object state to the given database.
     */
    void
    putToDatabaseSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBImplicitStaggeredHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBImplicitStaggeredHierarchyIntegrator(
        const IBImplicitStaggeredHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBImplicitStaggeredHierarchyIntegrator&
    operator=(
        const IBImplicitStaggeredHierarchyIntegrator& that);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void
    getFromRestart();

    double d_current_time, d_new_time;
    SAMRAI::tbox::Pointer<INSStaggeredStokesOperator> d_stokes_op;

    /*!
     * \brief Implementation of a modified Stokes operator that includes forcing
     * terms computed via the IB method.
     */
    class Operator
        : public IBTK::GeneralOperator
    {
    public:
        /*!
         * \brief Constructor.
         */
        Operator(
            const IBImplicitStaggeredHierarchyIntegrator* ib_solver);

        /*!
         * \brief Destructor.
         */
        ~Operator();

        /*!
         * \name General operator functionality.
         */
        //\{

        /*!
         * \brief Compute \f$y=F[x]\f$.
         */
        void
        apply(
            SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
            SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

        /*!
         * \brief Compute hierarchy dependent data required for computing y=F[x] and
         * z=F[x]+y.
         */
        void
        initializeOperatorState(
            const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
            const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out);

        /*!
         * \brief Remove all hierarchy dependent data allocated by
         * initializeOperatorState().
         */
        void
        deallocateOperatorState();

        //\}

        /*!
         * \name Logging functions.
         */
        //\{

        /*!
         * \brief Enable or disable logging.
         *
         * \param enabled logging state: true=on, false=off
         */
        void
        enableLogging(
            bool enabled=true);

        //\}

    private:
        /*!
         * \brief Default constructor.
         *
         * \note This constructor is not implemented and should not be used.
         */
        Operator();

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        Operator(
            const Operator& from);

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        Operator&
        operator=(
            const Operator& that);

        const IBImplicitStaggeredHierarchyIntegrator* const d_ib_solver;
    };

    /*!
     * \brief Implementation of the Jacobian of a modified Stokes operator that
     * includes forcing terms computed via the IB method.
     */
    class Jacobian
        : public IBTK::JacobianOperator
    {
    public:
        /*!
         * \brief Constructor.
         */
        Jacobian(
            const IBImplicitStaggeredHierarchyIntegrator* ib_solver);

        /*!
         * \brief Destructor.
         */
        ~Jacobian();

        /*!
         * \name General Jacobian functionality.
         */
        //\{

        /*!
         * \brief Compute hierarchy dependent data required for evaluating F'[x].
         *
         * \param x value where the Jacobian is to be evaluated
         */
        void
        formJacobian(
            SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x);

        /*!
         * \brief Return the vector where the Jacobian is evaluated.
         */
        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >
        getBaseVector() const;

        //\}

        /*!
         * \name Linear operator functionality.
         */
        //\{

        /*!
         * \brief Compute \f$y=F[x]\f$.
         */
        void
        apply(
            SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
            SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

        /*!
         * \brief Compute hierarchy dependent data required for computing y=F[x] and
         * z=F[x]+y.
         */
        void
        initializeOperatorState(
            const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
            const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out);

        /*!
         * \brief Remove all hierarchy dependent data allocated by
         * initializeOperatorState().
         */
        void
        deallocateOperatorState();

        //\}

        /*!
         * \name Logging functions.
         */
        //\{

        /*!
         * \brief Enable or disable logging.
         *
         * \param enabled logging state: true=on, false=off
         */
        void
        enableLogging(
            bool enabled=true);

        //\}

    private:
        /*!
         * \brief Default constructor.
         *
         * \note This constructor is not implemented and should not be used.
         */
        Jacobian();

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        Jacobian(
            const Jacobian& from);

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        Jacobian&
        operator=(
            const Jacobian& that);

        const IBImplicitStaggeredHierarchyIntegrator* const d_ib_solver;
        Mat d_J_mat;
        bool d_J_is_set;
        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_x_base;
    };
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBImplicitStaggeredHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBImplicitStaggeredHierarchyIntegrator
