// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
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

#ifndef included_IBAMR_INSStaggeredStochasticForcing
#define included_IBAMR_INSStaggeredStochasticForcing

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIArray.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIDatabase.h"
#include "SAMRAIEdgeVariable.h" // IWYU pragma: keep
#include "SAMRAIIntVector.h"
#include "SAMRAINodeVariable.h" // IWYU pragma: keep
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableContext.h"

#include <string>
#include <vector>

namespace IBAMR
{
class INSStaggeredHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class Variable;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI
/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredStochasticForcing provides an interface for
 * specifying a stochastic forcing term for a staggered-grid incompressible
 * Navier-Stokes solver.
 */
class INSStaggeredStochasticForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief This constructor creates Variable and VariableContext objects for
     * storing the stochastic stresses at the centers and nodes of the Cartesian
     * grid.
     */
    INSStaggeredStochasticForcing(std::string object_name,
                                  SAMRAIPointer<SAMRAIDatabase> input_db,
                                  const INSStaggeredHierarchyIntegrator* fluid_solver);

    /*!
     * \brief Empty destructor.
     */
    ~INSStaggeredStochasticForcing() = default;

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete INSStaggeredStochasticForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 SAMRAIPointer<SAMRAIVariable> var,
                                 SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = IBTK::invalid_level_number,
                                 const int finest_ln = IBTK::invalid_level_number) override;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void
    setDataOnPatch(const int data_idx,
                   SAMRAIPointer<SAMRAIVariable> var,
                   SAMRAIPointer<SAMRAIPatch> patch,
                   const double data_time,
                   const bool initial_time = false,
                   SAMRAIPointer<SAMRAIPatchLevel> patch_level = SAMRAIPointer<SAMRAIPatchLevel>(nullptr)) override;

    //\}

protected:
    /*!
     * The object name is used for error/warning reporting.
     */
    std::string d_object_name;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredStochasticForcing() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredStochasticForcing(const INSStaggeredStochasticForcing& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredStochasticForcing& operator=(const INSStaggeredStochasticForcing& that) = delete;

    /*!
     * Pointer to the fluid solver object that is using this stochastic force
     * generator.
     */
    const INSStaggeredHierarchyIntegrator* const d_fluid_solver;

    /*!
     * Type of stress tensor (correlated or uncorrelated).
     */
    StochasticStressTensorType d_stress_tensor_type = UNCORRELATED;

    /*!
     * Weighting data.
     */
    double d_std = std::numeric_limits<double>::quiet_NaN();
    int d_num_rand_vals = 0;
    std::vector<SAMRAIArray<double> > d_weights;

    /*!
     * Boundary condition scalings.
     */
    double d_velocity_bc_scaling, d_traction_bc_scaling = 0.0;

    /*!
     * VariableContext and Variable objects for storing the components of the
     * stochastic stresses.
     */
    SAMRAIPointer<SAMRAIVariableContext> d_context;
    SAMRAIPointer<SAMRAICellVariable<double> > d_W_cc_var;
    int d_W_cc_idx = IBTK::invalid_index;
    std::vector<int> d_W_cc_idxs;
#if (NDIM == 2)
    SAMRAIPointer<SAMRAINodeVariable<double> > d_W_nc_var;
    int d_W_nc_idx = IBTK::invalid_index;
    std::vector<int> d_W_nc_idxs;
#endif
#if (NDIM == 3)
    SAMRAIPointer<SAMRAIEdgeVariable<double> > d_W_ec_var;
    int d_W_ec_idx = IBTK::invalid_index;
    std::vector<int> d_W_ec_idxs;
#endif
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSStaggeredStochasticForcing
