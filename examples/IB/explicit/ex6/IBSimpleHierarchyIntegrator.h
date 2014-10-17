// Filename: IBSimpleHierarchyIntegrator.h
// Created on 04 Apr 2013 by Boyce Griffith
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

#ifndef included_IBSimpleHierarchyIntegrator
#define included_IBSimpleHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class IBSimpleHierarchyIntegrator is an implementation of a simple
 * first-order accurate, semi-implicit version of the immersed boundary method.
 */
class IBSimpleHierarchyIntegrator : public IBHierarchyIntegrator
{
public:
    /*!
     * The constructor for class IBSimpleHierarchyIntegrator sets some default
     * values and reads in configuration information from input and restart
     * databases.
     *
     * \warning This simple example class does not support restarting.
     */
    IBSimpleHierarchyIntegrator(const std::string& object_name,
                                Pointer<Database> input_db,
                                Pointer<IBMethod> ib_method_ops,
                                Pointer<INSHierarchyIntegrator> ins_hier_integrator);

    /*!
     * The destructor for class IBSimpleHierarchyIntegrator does
     * not do anything interesting.
     */
    ~IBSimpleHierarchyIntegrator();

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1);

    /*!
     * Initialize any variables, communications algorithms, solvers, or other
     * data structures required by this time integrator object.
     */
    void initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                       Pointer<GriddingAlgorithm<NDIM> > gridding_alg);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBSimpleHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBSimpleHierarchyIntegrator(const IBSimpleHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBSimpleHierarchyIntegrator& operator=(const IBSimpleHierarchyIntegrator& that);

    /*
     * Pointers to Lagrangian data objects.
     */
    Pointer<LData> d_X_current_data, d_X_new_data, d_U_data, d_F_data;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBSimpleHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBSimpleHierarchyIntegrator
