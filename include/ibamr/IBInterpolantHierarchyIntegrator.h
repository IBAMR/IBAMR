// Filename: IBInterpolantHierarchyIntegrator.h
// Created on 20 Nov 2018 by Amneet Bhalla
//
// Copyright (c) 2002-2018, Amneet Bhalla
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

#ifndef included_IBInterpolantHierarchyIntegrator
#define included_IBInterpolantHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include "ibamr/IBHierarchyIntegrator.h"

namespace IBAMR
{
class IBFEMethod;
class IBLevelSetMethod;
class IBInterpolantMethod;
class INSHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBInterpolantHierarchyIntegrator is an implementation of Brinkman
 * penalization immersed boundary method.
 */
class IBInterpolantHierarchyIntegrator : public IBAMR::IBHierarchyIntegrator
{
public:
    /*!
     * The constructor for class IBInterpolantHierarchyIntegrator sets some default
     * values and reads in configuration information from input and restart
     * databases.
     */
    IBInterpolantHierarchyIntegrator(const std::string& object_name,
                                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                     SAMRAI::tbox::Pointer<IBAMR::IBLevelSetMethod> ib_ls_method_ops,
                                     SAMRAI::tbox::Pointer<INSHierarchyIntegrator> ins_hier_integrator,
                                     bool register_for_restart = true);

    /*!
     * The destructor for class IBInterpolantHierarchyIntegrator does
     * not do anything interesting.
     */
    ~IBInterpolantHierarchyIntegrator();

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
    void initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

protected:
    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBInterpolantHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBInterpolantHierarchyIntegrator(const IBInterpolantHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBInterpolantHierarchyIntegrator& operator=(const IBInterpolantHierarchyIntegrator& that);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*
     * Pointers to IBMethod ops and INSHierarchyIntegrator.
     */
    SAMRAI::tbox::Pointer<IBLevelSetMethod> d_ib_ls_method_ops;
    SAMRAI::tbox::Pointer<IBInterpolantMethod> d_ib_interpolant_method_ops;
    SAMRAI::tbox::Pointer<INSHierarchyIntegrator> d_ins_hier_integrator;
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBInterpolantHierarchyIntegrator
