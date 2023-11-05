// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_LaserBeamForceFunction
#define included_IBAMR_LaserBeamForceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/PhaseChangeHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyMathOps.h"

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "PatchLevel.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"

#include <string>

// IBAMR INCLUDES
#include <ibamr/app_namespaces.h>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
template <int DIM>
class Patch;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class SideData;
template <int DIM, class TYPE>
class CellData;
} // namespace pdat
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LaserBeamForceFunction computes the source term due to Laser beam for energy equation.
 *
 * \note Presently, this class assumes that the indicator function is a cell centered
 * Heaviside variable that is maintained by the advection-diffusion integrator.
 *
 * Reference
 * Zaki Saptari Saldi, <A
 * HREF="https://repository.tudelft.nl/islandora/object/uuid%3A8401374b-9e9c-4d25-86b7-fc445ec73d27"> Marangoni driven
 * free surface flows in liquid weld pools</A>
 */
class LaserBeamForceFunction : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    LaserBeamForceFunction(const std::string& object_name,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           SAMRAI::tbox::Pointer<PhaseChangeHierarchyIntegrator> phase_change_solver,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > ls_var,
                           const double rho_liquid,
                           const double rho_gas,
                           const double cp_liquid,
                           const double cp_gas);

    /*!
     * \brief Destructor.
     */
    virtual ~LaserBeamForceFunction() = default;

    /*!
     * \name Methods to set the data.
     */
    //\{

    /*!
     * \note This concrete IBTK::CartGridFunction is time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy using the virtual function
     * setDataOnPatch().
     *
     * \see setDataOnPatch
     */
    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1) override;

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL)) override;

    /*!
     * \brief Function to compute heat influx.
     */
    using HeatInfluxPtr = void (*)(int F_idx,
                                   SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                   int cycle_num,
                                   double time,
                                   double current_time,
                                   double new_time,
                                   void* ctx);

    /*!
     * \brief Register function to compute heat influx.
     */
    void registerHeatInflux(HeatInfluxPtr callback, void* ctx);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LaserBeamForceFunction() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LaserBeamForceFunction(const LaserBeamForceFunction& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LaserBeamForceFunction& operator=(const LaserBeamForceFunction& that) = delete;

    /*!
     * Pointer to phase change solver.
     */
    SAMRAI::tbox::Pointer<PhaseChangeHierarchyIntegrator> d_phase_change_solver;

    /*!
     * Convert the level set variable to a smoothed heaviside function.
     */
    void convertToHeaviside(int H_idx,
                            int phi_idx,
                            int coarsest_ln,
                            int finest_ln,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * Mollify data.
     */
    void mollifyData(int H_scratch_idx,
                     int coarsest_ln,
                     int finest_ln,
                     double data_time,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                     SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> fill_op);

    /*!
     *  Variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_phi_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_grad_H_var;

    /*!
     * Scratch data.
     */
    int d_H_scratch_idx = IBTK::invalid_index, d_grad_H_scratch_idx = IBTK::invalid_index;

    /*!
     * Function to set heat influx from the main driver.
     */
    std::vector<HeatInfluxPtr> d_heat_influx;
    std::vector<void*> d_heat_influx_ctx;

    /*!
     * Time stepping type.
     */
    TimeSteppingType d_ts_type;

    /*!
     * Density of liquid and gas phases.
     */
    const double d_rho_liquid, d_rho_gas;

    /*!
     * Specific heat of liquid and gas phases.
     */
    const int d_cp_liquid, d_cp_gas;

    /*!
     * Get the stencil size for the kernel.
     */
    int getStencilSize(const std::string& kernel_fcn);

    /*!
     * Get the ghost cell width of scratch data.
     */
    int getMinimumGhostWidth(const std::string& kernel_fcn);

    std::string d_kernel_fcn;

    /*!
     * Number of interface cells.
     */
    double d_num_interface_cells;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_LaserBeamForceFunction
