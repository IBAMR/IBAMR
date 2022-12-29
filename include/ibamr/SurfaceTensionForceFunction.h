// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
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

#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "PatchLevel.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"

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
                                const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > level_set_var);

    /*!
     * \brief Destructor.
     */
    virtual ~SurfaceTensionForceFunction() = default;

    /*!
     * \brief Set the smoother (kernel function) to mollify the Heaviside function.
     */
    virtual void setSmoother(const std::string& kernel_fcn);

    /*!
     * \brief Set the constant surface tension coefficient.
     */
    virtual void setSurfaceTensionCoef(double sigma);

    /*!
     * \brief Set the number of interface cells m over which the surface tension force will be applied.
     * The surface tension will take effect in the band -m*h to m*h around the interface,
     * where h = (dx*dy)^(1/2) in 2D and h = (dx*dy*dz)^(1/3) in 3D.
     */
    virtual void setNumberOfInterfaceCells(double m);

    /*!
     * \brief Get the smoother (kernel function) to mollify the Heaviside function.
     */
    std::string getSmoother() const
    {
        return d_kernel_fcn;
    } // getSmoother

    /*!
     * \brief Get the constant surface tension coefficient.
     */
    double getSurfaceTensionCoef() const
    {
        return d_sigma;
    } // getSurfaceTensionCoef

    /*!
     * \brief Get the number of interface cells over which the surface tension force will be applied.
     */
    double getNumberOfInterfaceCells() const
    {
        return d_num_interface_cells;
    } // getNumberOfInterfaceCells

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
                                 int coarsest_ln = IBTK::invalid_level_number,
                                 int finest_ln = IBTK::invalid_level_number) override;

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

    //\}

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

    /*!
     * Convert the level set variable to a smoothed heaviside function.
     */
    void convertToHeaviside(int phi_idx,
                            int coarsest_ln,
                            int finest_ln,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * Mollify data.
     */
    void mollifyData(int phi_idx,
                     int coarsest_ln,
                     int finest_ln,
                     double data_time,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                     SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> fill_op);

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
     * Get the stencil size for the kernel.
     */
    int getStencilSize(const std::string& kernel_fcn);

    /*!
     * Get the ghost cell width of scratch data.
     */
    int getMinimumGhostWidth(const std::string& kernel_fcn);

    const AdvDiffHierarchyIntegrator* const d_adv_diff_solver;
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_ls_var;
    TimeSteppingType d_ts_type;
    int d_C_idx, d_phi_idx;
    std::string d_kernel_fcn;
    double d_sigma, d_num_interface_cells;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_SurfaceTensionForceFunction
