// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2014 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

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
