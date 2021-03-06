// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
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

#ifndef included_IBHydrodynamicSurfaceForceEvaluator
#define included_IBHydrodynamicSurfaceForceEvaluator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/INSHierarchyIntegrator.h"

#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CellVariable.h"
#include "RobinBcCoefStrategy.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "Eigen/Core"
#include "Eigen/Geometry"
IBTK_ENABLE_EXTRA_WARNINGS

#include <iosfwd>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////
namespace IBAMR
{
class AdvDiffHierarchyIntegrator;
class INSHierarchyIntegrator;
} // namespace IBAMR
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
class Database;
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
class IBHydrodynamicSurfaceForceEvaluator : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IBHydrodynamicSurfaceForceEvaluator(std::string object_name,
                                        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                                        SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                                        SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> fluid_solver,
                                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db = NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBHydrodynamicSurfaceForceEvaluator();

    /*!
     * \brief Compute the hydrodynamic force via surface integration after all the hierarchy
     * variables have been integrated.
     *
     * \note The users should call this routine if the hydrodynamic forces need to be computed just
     * as a postprocessing step.
     *
     * @note This function uses 3D vectors even in 2D so that the torque makes sense.
     */
    virtual void computeHydrodynamicForceTorque(IBTK::Vector3d& pressure_force,
                                                IBTK::Vector3d& viscous_force,
                                                IBTK::Vector3d& pressure_torque,
                                                IBTK::Vector3d& viscous_torque,
                                                const IBTK::Vector3d& X0);

    /*!
     * \brief Compute the hydrodynamic force via surface integration while the integration of
     * variables is happening during a given timestep.
     *
     * @note This function uses 3D vectors even in 2D so that the torque makes sense.
     */
    virtual void computeHydrodynamicForceTorque(IBTK::Vector3d& pressure_force,
                                                IBTK::Vector3d& viscous_force,
                                                IBTK::Vector3d& pressure_torque,
                                                IBTK::Vector3d& viscous_torque,
                                                const IBTK::Vector3d& X0,
                                                double time,
                                                double current_time,
                                                double new_time);
    /*!
     * \brief Set the surface contour value <em> s <\em> of the level set field, such that \f$ \phi(s) = \partial \Omega
     * \f$.
     */
    virtual void setSurfaceContourLevel(double s = 0);

    /*!
     * \brief Indicate if the force and torque results need to be written on a file.
     */
    void writeToFile(bool write_to_file = true);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBHydrodynamicSurfaceForceEvaluator(const IBHydrodynamicSurfaceForceEvaluator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBHydrodynamicSurfaceForceEvaluator& operator=(const IBHydrodynamicSurfaceForceEvaluator& that) = delete;

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Fill required patch data and ghost cells.
     */
    void fillPatchData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                       double fill_time,
                       bool use_current_ctx,
                       bool use_new_ctx);

    /*!
     * \brief Object name.
     */
    std::string d_object_name;

    /*!
     * \brief Level set variable for the immersed body
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_solid_var;

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
    int d_ls_solid_idx = IBTK::invalid_index;

    /*!
     * \brief Fluid velocity and pressure patch data indices.
     */
    int d_u_idx = IBTK::invalid_index, d_p_idx = IBTK::invalid_index;

    /*!
     * \brief Viscosity patch data index.
     */
    int d_mu_idx = IBTK::invalid_index;

    /*!
     * \brief Viscosity for the constant mu case.
     */
    double d_mu = std::numeric_limits<double>::quiet_NaN();

    /*!
     * \brief Whether or not viscosity is constant.
     */
    bool d_mu_is_const;

    /*!
     * \brief The contour level that describes the surface of the solid object.
     */
    double d_surface_contour_value = 0.0;

    /*!
     * \brief Whether to write results on a text file.
     */
    bool d_write_to_file = false;

    /*!
     * \brief File streams associated for the output of hydrodynamic force.
     *
     * \note Columns 1-3 represent sum of -p.n dA. Columns 4-6 represent sum of n.(grad U + grad U^T) dA.
     *
     */
    std::unique_ptr<std::ofstream> d_hydro_force_stream = nullptr;

    /*!
     * \brief File streams associated for the output of hydrodynamic torque.
     *
     * \note Columns 1-3 represent sum of r X -p.n dA. Columns 4-6 represent sum of r x n.(grad U + grad U^T) dA.
     *
     */
    std::unique_ptr<std::ofstream> d_hydro_torque_stream = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBHydrodynamicSurfaceForceEvaluator
