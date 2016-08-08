// Filename: ForceProjector.h
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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
// POSSIBILITY OF SUCH DAMAGE

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_ForceProjector
#define included_ForceProjector

///////////////////////////// INCLUDES ///////////////////////////////////

#include <string>

#include <PatchHierarchy.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

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

    //////////////// PRIVATE /////////////////////////////

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

} // IBTK

#endif // #ifndef included_ForceProjector
