// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2026 by the IBAMR developers
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

#ifndef included_IBAMR_SurfaceTensionForceFunction
#define included_IBAMR_SurfaceTensionForceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/ibamr_enums.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

#include <tbox/Array.h>
#include <tbox/Pointer.h>

#include <CartesianGridGeometry.h>
#include <IntVector.h>
#include <PatchLevel.h>

#include <string>



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
template <int DIM>
class PatchHierarchy;
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
 * \brief Class SurfaceTensionForceFunction provides surface tension forcing
 * using the continuum surface tension force model of Brackbill, Kothe, and Zemach.
 *
 * \note Presently, this class assumes that the indicator function is a cell centered
 * level-set variable that is maintained by the advection-diffusion integrator. In general,
 * the indicator variable can either be a level set function, a volume fraction function,
 * or a phase field function.
 *
 * Reference
 * Brackbill et. al, <A HREF="https://www.sciencedirect.com/science/article/pii/002199919290240Y">
 * A continuum method for modeling surface tension</A>
 */
class SurfaceTensionForceFunction : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    SurfaceTensionForceFunction(const std::string& object_name,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                const AdvDiffHierarchyIntegrator* adv_diff_solver,
                                const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> level_set_var);

    /*!
     * \brief Destructor.
     */
     ~SurfaceTensionForceFunction() override = default;



    /*!
     * \brief Set the smoother (kernel function) to mollify the Heaviside function.
     */
    //  void setSmoother(const std::string& kernel_fcn);

    /*!
     * \brief Set the constant surface tension coefficient.
     */
    virtual void setSurfaceTensionCoef(double sigma) = 0; // make it virtual


    /*!
     * \brief Get the constant surface tension coefficient.
     */
    virtual double getSurfaceTensionCoef() const = 0; //make it virtual

    /*!
     * \note This concrete IBTK::CartGridFunction is time-dependent.
     */
    bool isTimeDependent() const override = 0;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy using the virtual function
     * setDataOnPatch().
     *
     * \see setDataOnPatch
     */
    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = IBTK::invalid_level_number,
                                 int finest_ln = IBTK::invalid_level_number) override= 0;

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(nullptr)) override= 0;

    /*!
     * \brief Function to Mask surface tension force to act only on the liquid-gas interface.
     */
    using MaskSurfaceTensionForcePtr = void (*)(int F_idx,
                                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                int cycle_num,
                                                double time,
                                                double current_time,
                                                double new_time,
                                                void* ctx);

    /*!
     * \brief Register function to limit the surface tension force.
     */
    virtual void registerSurfaceTensionForceMasking(MaskSurfaceTensionForcePtr callback, void* ctx) = 0;

    /*!
     * \brief Function to compute the variable surface tension coefficient.
     */
    using ComputeSurfaceTensionCoefficientPtr = void (*)(int F_idx,
                                                         SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                         int cycle_num,
                                                         double time,
                                                         double current_time,
                                                         double new_time,
                                                         void* ctx);

    /*!
     * \brief Register function to compute the variable surface tension coefficient.
     */
    virtual void registerSurfaceTensionCoefficientFunction(ComputeSurfaceTensionCoefficientPtr callback, void* ctx) = 0;

    //\}
protected:
   
        // common to both the derived classes
    const AdvDiffHierarchyIntegrator* const d_adv_diff_solver;

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> d_ls_var;
    

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SurfaceTensionForceFunction() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SurfaceTensionForceFunction(const SurfaceTensionForceFunction& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SurfaceTensionForceFunction& operator=(const SurfaceTensionForceFunction& that) = delete;

};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_SurfaceTensionForceFunction
