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

#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
// #include "ibamr/AdvDiffHierarchyIntegrator.h"
// #include "ibamr/IEPSemiImplicitHierarchyIntegrator.h"

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "PatchLevel.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"

#include <string>

// IBAMR INCLUDES
#include <ibamr/app_namespaces.h>

namespace IBAMR
{
class AdvDiffHierarchyIntegrator;
class IEPSemiImplicitHierarchyIntegrator;
} // namespace IBAMR
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
                           // AdvDiffHierarchyIntegrator* adv_diff_solver,
                           // SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                           SAMRAI::tbox::Pointer<IEPSemiImplicitHierarchyIntegrator> iep_solver,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > H_var,
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

    // AdvDiffHierarchyIntegrator* d_adv_diff_solver;
    // SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<IEPSemiImplicitHierarchyIntegrator> d_iep_solver;
    // IEPSemiImplicitHierarchyIntegrator* d_iep_solver;
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_H_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_grad_H_var;
    TimeSteppingType d_ts_type;
    int d_H_scratch_idx, d_grad_H_scratch_idx;

    std::vector<HeatInfluxPtr> d_heat_influx;
    std::vector<void*> d_heat_influx_ctx;

    const double d_rho_liquid, d_rho_gas, d_cp_liquid, d_cp_gas;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_LaserBeamForceFunction
