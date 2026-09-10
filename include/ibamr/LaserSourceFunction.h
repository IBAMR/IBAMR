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

#ifndef included_IBAMR_LaserSourceFunction
#define included_IBAMR_LaserSourceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/PhaseChangeHierarchyIntegrator.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyMathOps.h>

#include <tbox/Array.h>
#include <tbox/Pointer.h>

#include <CartesianGridGeometry.h>
#include <IntVector.h>
#include <PatchLevel.h>

#include <string>

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
 * \brief Apply an interfacial heat flux as a cell-centered energy source.
 *
 * Supply a cell-centered level set registered with the phase-change integrator;
 * its nonnegative side denotes the heated material. The heat-flux callback scales the interfacial weights in place. The
 * solver and variable are retained; the callback context is borrowed and must remain valid during evaluation.
 *
 * Call setDataOnPatchHierarchy() to prepare the level-set-derived interfacial
 * weights before patch evaluation. Initial-time evaluation supplies zero weights to the callback.
 *
 * Reference
 * Thirumalaisamy and Bhalla, <A
 * HREF="https://arxiv.org/abs/2407.05588"> A consistent, volume preserving, and adaptive mesh
 * refinement-based framework for modeling non-isothermal gas-liquid-solid flows with phase change</A>
 */
class LaserSourceFunction : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    LaserSourceFunction(const std::string& object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        SAMRAI::tbox::Pointer<PhaseChangeHierarchyIntegrator> phase_change_solver,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> level_set_var);

    /*!
     * \brief Destructor.
     */
    ~LaserSourceFunction() override = default;

    /*!
     * \name Methods to set the data.
     */
    //\{

    /*!
     * \note This concrete IBTK::CartGridFunction is time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Form interfacial weights and apply the registered heat flux.
     */
    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1) override;

    /*!
     * \brief Evaluate the magnitude of the prepared Heaviside gradient on a patch.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(nullptr)) override;

    /*!
     * \brief Callback function to compute the imposed heat flux.
     */
    using HeatFluxPtr = void (*)(int F_idx,
                                 SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                 int cycle_num,
                                 double time,
                                 double current_time,
                                 double new_time,
                                 void* ctx);

    /*!
     * \brief Register function to compute the imposed heat flux.
     */
    void registerHeatFlux(HeatFluxPtr callback, void* ctx);

private:
    LaserSourceFunction() = delete;

    LaserSourceFunction(const LaserSourceFunction& from) = delete;

    LaserSourceFunction& operator=(const LaserSourceFunction& that) = delete;

    /*!
     * Convert the level set variable to a smoothed Heaviside or discontinuous Heaviside function.
     */
    void convertToHeaviside(int H_idx,
                            int phi_idx,
                            int coarsest_ln,
                            int finest_ln,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy);

    /*!
     * Mollify data.
     */
    void mollifyData(int H_scratch_idx,
                     int coarsest_ln,
                     int finest_ln,
                     double data_time,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                     SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> fill_op);

    /*!
     * Pointer to phase change solver.
     */
    SAMRAI::tbox::Pointer<PhaseChangeHierarchyIntegrator> d_phase_change_solver;

    /*!
     *  Variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> d_phi_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_grad_H_var;

    /*!
     * Scratch data.
     */
    int d_H_scratch_idx = IBTK::invalid_index, d_grad_H_scratch_idx = IBTK::invalid_index;

    /*!
     * Call back function and the context to find the heat flux.
     */
    HeatFluxPtr d_heat_flux = nullptr;
    void* d_heat_flux_ctx = nullptr;

    /*!
     * Time stepping type.
     */
    TimeSteppingType d_ts_type;

    std::string d_kernel_fcn;

    /*!
     * Number of interface cells.
     */
    double d_num_interface_cells = std::numeric_limits<double>::signaling_NaN();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_LaserSourceFunction
