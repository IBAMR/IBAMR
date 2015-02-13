#ifndef included_NonbondedInteractions
#define included_NonbondedInteractions

// IBAMR INCLUDES
#include "tbox/Array.h"
#include "ArrayData.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LNodeSetData.h"
#include "ibamr/IBSpringForceSpec.h"
#include <ibtk/IndexUtilities.h>
#include <ibtk/AppInitializer.h>
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/app_namespaces.h>
#include <ibamr/IBLagrangianForceStrategy.h>


// IBTK THIRD-PARTY INCLUDES
#include <muParser.h>


// Nonbonded Force Function Pointer.
// Takes vector between
// the points, D = q_i - q_j and sets values of
// out_force, the force that q_i experiences.
// q_j will experience force -out_force.
//
// parameters are passed in the double* params.
typedef void (*NonBddForceFcnPtr)(double* D,
                                  const SAMRAI::tbox::Array<double> params,
                                  double* out_force);

/////////////////////////////// NonbondedForceEvaluator ////////////////////////

class NonbondedForceEvaluator : public IBLagrangianForceStrategy
{

  public:
    
    // constructor
    NonbondedForceEvaluator(Pointer<Database> input_db,
                            Pointer<CartesianGridGeometry<NDIM> > grid_geometry);
    
    // function to evaluate forces
    void
    EvaluateForces(int mstr_petsc_idx,
                   int search_petsc_idx,
                   Pointer<LData> X_data,
                   vector<int> cell_offset,
                   Pointer<LData> F_data);

    // Overwrite the computeLagrangianForce function.
    void
    computeLagrangianForce(Pointer<LData> F_data,
                           Pointer<LData> X_data,
                           Pointer<LData> U_data,
                           const Pointer<PatchHierarchy<NDIM> > hierarchy,
                           const int level_number,
                           const double data_time,
                           LDataManager* const l_data_manager);

    // register the inline function used
    void
    registerForceFcnPtr(NonBddForceFcnPtr force_fcn_ptr);

  private:
    
    // default constructor, not implemented
    NonbondedForceEvaluator();

    // copy constructor, not implemented
    NonbondedForceEvaluator(const NonbondedForceEvaluator& from);

    // assignment operator, not implemented
    NonbondedForceEvaluator&
    operator=(const NonbondedForceEvaluator& that);
    

    // which type of force do we use?
    int d_force_type;

    // interaction radius
    double d_interaction_radius;
    
    // regrid_alpha, for computing buffer to add to interactions
    double d_regrid_alpha;

    // parameters for force function.
    SAMRAI::tbox::Array<double> d_parameters;
    
    // grid geometry
    Pointer<CartesianGridGeometry<NDIM> > d_grid_geometry;

    // spring force function pointer.  Used to evaluate the force between particles
    //  TODO: Add species, make this a map from species1 x species2 -> Force Function Pointer
    NonBddForceFcnPtr d_force_fcn_ptr;

};


#endif
