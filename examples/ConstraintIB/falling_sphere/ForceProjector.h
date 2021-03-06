// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_ForceProjector
#define included_ForceProjector

///////////////////////////// INCLUDES ///////////////////////////////////

#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <PatchHierarchy.h>
#include <SideVariable.h>
#include <VariableContext.h>

#include <string>

namespace IBTK
{
/*!
 * Pre processing call back function to be hooked into IBAMR::INSStaggeredHierachyIntegrator class.
 *
 * \param current_time is the time at t_n.
 * \param new_time is the time at t_n+1 = t_n + dt.
 * \param cycle_num is the cycle of predictor-corrector scheme.
 * \param ctx is the pointer to IBTK::ForceProjector class object.
 */

void
callForceProjectorCallBackFunction(const double current_time, const double new_time, const int cycle_num, void* ctx);

/*!
 * \brief Class ForceProjector is a utility class which projects force from
 * Lagrangian points onto the background mesh.
 *
 */
class ForceProjector
{
public:
    /*!
     * \brief Constructor.
     */
    ForceProjector(const std::string& object_name,
                   IBTK::LDataManager* lag_data_manager,
                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                   const std::string solver_type = "STAGGERED");

    /*!
     * \brief Destructor.
     */
    ~ForceProjector();

    /*!
     * \brief Register the name of Lagrangian quantities to be used to calculate forces on Lagrangian
     * points.
     */
    void registerLagrangianQuantityName(const std::string& lag_quantity_name);

    void registerLagrangianQuantitiesName(const std::vector<std::string>& lag_quantities_name);

    /*!
     * \brief Register volume associated with each material point.
     */
    void associateVolumeElement(const double vol_lag_pt);

    /*!
     * \brief Calculate forcing on Lagrangian points.
     */
    void calculateLagrangianBodyForce(const double new_time, const double current_time);

    /*!
     * \brief Spread the Lagrangian forcing on the background mesh.
     */
    void calculateEulerianBodyForce(const double new_time, const double current_time);

    /*!
     *  \brief Get the patch index associated with Eulerian force.
     */
    inline int getEulerianForcePatchDataIndex() const
    {
        return d_body_force_idx;
    } // getEulerianForcePatchDataIndex

private:
    /*!
     * \brief Default constructor is not implemented and should not be used.
     */
    ForceProjector();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    ForceProjector& operator=(const ForceProjector& that);

    /*!
     * \brief Default copy constructor is not implemented and should not be used.
     */
    ForceProjector(const ForceProjector& from);

    /*!
     * \brief Get the values from input_db.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to LDataManager.
     */
    IBTK::LDataManager* d_lag_data_manager;

    /*!
     * Pointer to Patch Hierarchy.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_patch_hierarchy;

    /*!
     * Fluid solver type.
     */
    std::string d_solver_type;

    /*!
     * Pointer to Lagrangian force data.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_lag_force;

    /*!
     *  Variables and variable context associated with calculating Eulerian force.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_body_force_var;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_body_force_context;
    int d_body_force_idx;

    /*!
     * Name of Lagrangian quantities to be used in calculating forces.
     */
    std::vector<std::string> d_lag_quantities_name;

    /*!
     * Volume associated with each element.
     */
    double d_vol_lag_pt;

    /*!
     * Densities of fluid and body.
     */
    double d_rho_fluid, d_rho_body;

    /*!
     * Gravitational force constants.
     */
    SAMRAI::tbox::Array<double> d_grav_const;

}; // ForceProjector

} // namespace IBTK

#endif // #ifndef included_ForceProjector
