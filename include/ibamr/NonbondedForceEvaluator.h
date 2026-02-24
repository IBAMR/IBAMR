// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2025 by the IBAMR developers
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

#ifndef included_IBAMR_NonbondedForceEvaluator
#define included_IBAMR_NonbondedForceEvaluator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/IBLagrangianForceStrategy.h>

#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIArray.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAIDatabase.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPointer.h>

#include <vector>

namespace IBTK
{
class LData;
class LDataManager;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

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
    using NonBddForceFcnPtr = void (*)(double* D, const SAMRAIArray<double> params, double* out_force);

    // Class constructor.
    NonbondedForceEvaluator(SAMRAIPointer<SAMRAIDatabase> input_db,
                            SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry);

    // Function to evaluate forces.
    void evaluateForces(int mstr_petsc_idx,
                        int search_petsc_idx,
                        SAMRAIPointer<IBTK::LData> X_data,
                        std::vector<int> cell_offset,
                        SAMRAIPointer<IBTK::LData> F_data);

    // Implementation of computeLagrangianForce.
    void computeLagrangianForce(SAMRAIPointer<IBTK::LData> F_data,
                                SAMRAIPointer<IBTK::LData> X_data,
                                SAMRAIPointer<IBTK::LData> U_data,
                                const SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
                                const int level_number,
                                const double data_time,
                                IBTK::LDataManager* const l_data_manager) override;

    // Register the force function used
    void registerForceFcnPtr(NonBddForceFcnPtr force_fcn_ptr);

private:
    // Default constructor, not implemented.
    NonbondedForceEvaluator() = delete;

    // Copy constructor, not implemented.
    NonbondedForceEvaluator(const NonbondedForceEvaluator& from) = delete;

    // Assignment operator, not implemented.
    NonbondedForceEvaluator& operator=(const NonbondedForceEvaluator& that) = delete;

    // interaction radius:
    double d_interaction_radius;

    // regrid_alpha, for computing buffer to add to interactions:
    double d_regrid_alpha;

    // parameters for force function:
    SAMRAIArray<double> d_parameters;

    // grid geometry
    SAMRAIPointer<SAMRAICartesianGridGeometry> d_grid_geometry;

    // spring force function pointer, to evaluate the force between particles:
    // TODO: Add species, make this a map from species1 x species2 -> Force Function Pointer
    NonBddForceFcnPtr d_force_fcn_ptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_NonbondedForceEvaluator
