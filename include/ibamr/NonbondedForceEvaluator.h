// Filename: NonbondexForceEvaluator.h
// Created by Steven Delong
//
// Copyright (c) 2002-2014, Steven Delong, Boyce Griffith
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

#ifndef included_NonbondedForceEvaluator
#define included_NonbondedForceEvaluator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "tbox/Array.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibamr/IBLagrangianForceStrategy.h"
#include "muParser.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
class NonbondedForceEvaluator : public IBLagrangianForceStrategy
{
public:
    // Nonbonded Force Function Pointer.
    // Takes vector between
    // the points, D = q_i - q_j and sets values of
    // out_force, the force that q_i experiences.
    // q_j will experience force -out_force.
    //
    // parameters are passed in the double* params.
    typedef void (*NonBddForceFcnPtr)(double* D, const SAMRAI::tbox::Array<double> params, double* out_force);

    // Class constructor.
    NonbondedForceEvaluator(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry);

    // Function to evaluate forces.
    void evaluateForces(int mstr_petsc_idx,
                        int search_petsc_idx,
                        SAMRAI::tbox::Pointer<IBTK::LData> X_data,
                        std::vector<int> cell_offset,
                        SAMRAI::tbox::Pointer<IBTK::LData> F_data);

    // Implementation of computeLagrangianForce.
    void computeLagrangianForce(SAMRAI::tbox::Pointer<IBTK::LData> F_data,
                                SAMRAI::tbox::Pointer<IBTK::LData> X_data,
                                SAMRAI::tbox::Pointer<IBTK::LData> U_data,
                                const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                const int level_number,
                                const double data_time,
                                IBTK::LDataManager* const l_data_manager);

    // Register the force function used
    void registerForceFcnPtr(NonBddForceFcnPtr force_fcn_ptr);

private:
    // Default constructor, not implemented.
    NonbondedForceEvaluator();

    // Copy constructor, not implemented.
    NonbondedForceEvaluator(const NonbondedForceEvaluator& from);

    // Assignment operator, not implemented.
    NonbondedForceEvaluator& operator=(const NonbondedForceEvaluator& that);

    // type of force to use:
    int d_force_type;

    // interaction radius:
    double d_interaction_radius;

    // regrid_alpha, for computing buffer to add to interactions:
    double d_regrid_alpha;

    // parameters for force function:
    SAMRAI::tbox::Array<double> d_parameters;

    // grid geometry
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geometry;

    // spring force function pointer, to evaluate the force between particles:
    // TODO: Add species, make this a map from species1 x species2 -> Force Function Pointer
    NonBddForceFcnPtr d_force_fcn_ptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_NonbondedForceEvaluator
