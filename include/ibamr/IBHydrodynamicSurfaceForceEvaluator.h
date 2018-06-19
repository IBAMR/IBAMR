// Filename: IBHydrodynamicSurfaceForceEvaluator.h
// Created on 11 May 2018 by Nishant Nangia
//
// Copyright (c) 2002-2018, Nishant Nangia and Amneet Bhalla
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

#ifndef included_IBHydrodynamicSurfaceForceEvaluator
#define included_IBHydrodynamicSurfaceForceEvaluator

#include <map>
#include <vector>

#include "Box.h"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "RobinBcCoefStrategy.h"
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSHierarchyIntegrator.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/ibtk_utilities.h>

namespace IBTK
{
class LData;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace tbox
{
template <class TYPE>
class Pointer;
} // namespace tbox
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBHydrodynamicSurfaceForceEvaluator computes hydrodynamic force
 * on immersed bodies. The class uses a level set indicator for the solid body to
 * add up surface forces in the cells adjacent to the immersed structure
 *
 * References
 * Patel, J,K. and Natarajan, G., <A HREF="https://www.sciencedirect.com/science/article/pii/S0021999118300342">Diffuse
 * interface
 * immersed boundary method for multi-fluid flows with arbitrarily moving rigid bodies.</A>
 *
 * \note  This class will work with both the constant and variable viscosity
 */
class IBHydrodynamicSurfaceForceEvaluator
{
public:
    /*!
     * \brief Default constructor.
     */
    IBHydrodynamicSurfaceForceEvaluator(const std::string& object_name,
                                        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                                        SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                                        SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> fluid_solver,
                                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db = NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBHydrodynamicSurfaceForceEvaluator();

    /*!
     * \brief Compute the hydrodynamic force via surface integration.
     */
    virtual void computeHydrodynamicForce();

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBHydrodynamicSurfaceForceEvaluator(const IBHydrodynamicSurfaceForceEvaluator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBHydrodynamicSurfaceForceEvaluator& operator=(const IBHydrodynamicSurfaceForceEvaluator& that);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Fill required patch data and ghost cells.
     */
    void fillPatchData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                       const double fill_time);

    /*!
     * \brief Object name.
     */
    std::string d_object_name;

    /*!
     * \brief Level set variable for the immersed body
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >  d_ls_solid_var;

    /*!
     * \brief Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * \brief Pointer to incompressible Navier-Stokes solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> d_fluid_solver;

    /*!
     * \brief Level set patch data index.
     */
    int d_ls_solid_idx;

    /*!
     * \brief Fluid velocity and pressure patch data indices.
     */
    int d_u_idx, d_p_idx;

    /*!
     * \brief Viscosity patch data index.
     */
    int d_mu_idx;

    /*!
     * \brief Viscosity for the constant mu case.
     */
    double d_mu;

    /*!
     * \brief Whether or not viscosity is constant.
     */
    bool d_mu_is_const;

    /*!
     * \brief The contour level that describes the surface of the solid object.
     */
    double d_surface_contour_value;

    /*!
     * \brief File streams associated for the output.
     *
     * \note Columns 1-3 represent sum of -p.n dA. Columns 4-6 represent sum of n.(grad U + grad U^T) dA.
     *
     */
    std::ofstream* d_hydro_force_stream;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBHydrodynamicSurfaceForceEvaluator
