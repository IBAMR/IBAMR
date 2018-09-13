// Filename: SurfaceTensionForceFunction.h
// Created on 16 Dec 2017 by Amneet Bhalla
//
// Copyright (c) 2002-2017, Amneet Bhalla
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBAMR_SurfaceTensionForceFunction
#define included_IBAMR_SurfaceTensionForceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "PatchLevel.h"
#include "boost/array.hpp"
#include "ibamr/ibamr_enums.h"
#include "ibtk/CartGridFunction.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"

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
    virtual ~SurfaceTensionForceFunction();

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
    bool isTimeDependent() const;

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
                                 int finest_ln = -1);

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SurfaceTensionForceFunction();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SurfaceTensionForceFunction(const SurfaceTensionForceFunction& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SurfaceTensionForceFunction& operator=(const SurfaceTensionForceFunction& that);

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
