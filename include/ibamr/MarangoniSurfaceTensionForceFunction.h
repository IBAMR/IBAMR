// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
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

#ifndef included_IBAMR_MarangoniSurfaceTensionForceFunction
#define included_IBAMR_MarangoniSurfaceTensionForceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/SurfaceTensionForceFunction.h"

namespace IBAMR
{
class AdvDiffHierarchyIntegrator;
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
 * \brief Class MarangoniSurfaceTensionForceFunction provides Marangoni forcing due to temperature
 * variations. This class computes the forcing term
 * \f$ F = \frac{\mathrm{d} \sigma}{\mathrm{d} T}(\nabla T |\nabla C| - (\nabla T \cdot \nabla \phi) \nabla C)\f$
 *  and add it to the surface tension forcing.
 *
 * This class uses the callback function to compute the variable Marangoni coefficient
 * as a function of Temperature.
 *
 * \note Presently, this class assumes that the indicator function is a cell centered
 * level-set variable that is maintained by the advection-diffusion integrator. In general,
 * the indicator variable can either be a level set function, a volume fraction function,
 * or a phase field function.
 */
class MarangoniSurfaceTensionForceFunction : public SurfaceTensionForceFunction
{
public:
    /*!
     * \brief Constructor.
     */
    MarangoniSurfaceTensionForceFunction(const std::string& object_name,
                                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                         AdvDiffHierarchyIntegrator* adv_diff_solver,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > level_set_var,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > T_var,
                                         SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& T_bc_coef);

    /*!
     * \brief Destructor.
     */
    virtual ~MarangoniSurfaceTensionForceFunction() = default;

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
     * \brief Callback function to compute the marangoni coefficient as a function of temperature and
     * multiply it with the F_data as F_data = marangoni_coef*F_data.
     */
    using ComputeMarangoniCoefPtr = void (*)(int F_idx,
                                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                             int cycle_num,
                                             double time,
                                             double current_time,
                                             double new_time,
                                             void* ctx);

    /*!
     * \brief Register callback function to compute the variable marangoni coefficient.
     */
    void registerMarangoniCoefficientFunction(ComputeMarangoniCoefPtr callback, void* ctx);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    MarangoniSurfaceTensionForceFunction() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    MarangoniSurfaceTensionForceFunction(const SurfaceTensionForceFunction& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    MarangoniSurfaceTensionForceFunction& operator=(const SurfaceTensionForceFunction& that) = delete;

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatchCell(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > F_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const double data_time,
                            const bool initial_time,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level);

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatchSide(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > F_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const double data_time,
                            const bool initial_time,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level);

    /*!
     * Temperature variable and its patch data index.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_T_var;
    int d_T_idx = IBTK::invalid_index, d_F_cloned_idx = IBTK::invalid_index;

    /*!
     * Marangoni coefficient.
     */
    double d_marangoni_coefficient = 0.0;

    /*!
     * Boundary condition object for temperature.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_T_bc_coef = nullptr;

    /*!
     * Call back function and the context to find marangoni coefficient as a function of temperature.
     */
    ComputeMarangoniCoefPtr d_compute_marangoni_coef = nullptr;
    void* d_compute_marangoni_coef_ctx = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_SurfaceTensionForceFunction
