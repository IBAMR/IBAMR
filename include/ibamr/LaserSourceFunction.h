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
 * \brief Class LaserSourceFunction can be used to account the laser beam effects as a source term into the energy
 * equation.
 *
 * \note Presently, this class assumes that the indicator function is a cell centered
 * Heaviside variable that is maintained by the advection-diffusion integrator.
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
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > level_set_var);

    /*!
     * \brief Destructor.
     */
    virtual ~LaserSourceFunction() = default;

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
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LaserSourceFunction() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LaserSourceFunction(const LaserSourceFunction& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LaserSourceFunction& operator=(const LaserSourceFunction& that) = delete;

    /*!
     * Convert the level set variable to a smoothed Heaviside or discontinuous Heaviside function.
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
     * Pointer to phase change solver.
     */
    SAMRAI::tbox::Pointer<PhaseChangeHierarchyIntegrator> d_phase_change_solver;

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
     * Call back function and the context to find the heat flux.
     */
    HeatFluxPtr d_heat_flux = nullptr;
    void* d_heat_flux_ctx = nullptr;

    /*!
     * Time stepping type.
     */
    TimeSteppingType d_ts_type;

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
    double d_num_interface_cells = std::numeric_limits<double>::signaling_NaN();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_LaserSourceFunction
