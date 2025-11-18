// Filename: CustomForceFunction.h
// Created by: Vinod
//
// Custom force function for applying forces to IB structures

#ifndef included_CustomForceFunction
#define included_CustomForceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/IBLagrangianForceStrategy.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <string>

namespace IBTK
{
class LDataManager;
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
class Database;
} // namespace tbox
} // namespace SAMRAI

/*!
 * \brief Class CustomForceFunction provides a custom implementation of
 * IBLagrangianForceStrategy for applying time-dependent forces to IB structures.
 *
 * This class demonstrates how to:
 * - Implement custom forcing strategies
 * - Apply position and time-dependent forces
 * - Access Lagrangian data in IBAMR
 *
 * Example usage patterns:
 * - Apply external forces (e.g., gravity, actuation)
 * - Implement feedback control
 * - Apply spatially varying forces
 */
class CustomForceFunction : public IBAMR::IBLagrangianForceStrategy
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name Name of the object
     * \param input_db Input database (optional)
     */
    CustomForceFunction(const std::string& object_name,
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr);

    /*!
     * \brief Destructor.
     */
    virtual ~CustomForceFunction();

    /*!
     * \brief Set the force magnitude.
     */
    void setForceMagnitude(double force_mag);

    /*!
     * \brief Set the force direction.
     */
    void setForceDirection(double fx, double fy);

    /*!
     * \brief Compute the force at the specified time on the Lagrangian structure.
     *
     * \param F_data Force data
     * \param X_data Position data
     * \param U_data Velocity data (optional, may be nullptr)
     * \param hierarchy Patch hierarchy
     * \param level_number Level number
     * \param data_time Current time
     * \param lag_manager Lagrangian data manager
     */
    virtual void setLagrangianForce(SAMRAI::tbox::Pointer<IBTK::LData> F_data,
                                   SAMRAI::tbox::Pointer<IBTK::LData> X_data,
                                   SAMRAI::tbox::Pointer<IBTK::LData> U_data,
                                   const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                   const int level_number,
                                   const double data_time,
                                   IBTK::LDataManager* const lag_manager);

private:
    /*!
     * \brief Copy constructor (not implemented).
     */
    CustomForceFunction(const CustomForceFunction& from) = delete;

    /*!
     * \brief Assignment operator (not implemented).
     */
    CustomForceFunction& operator=(const CustomForceFunction& that) = delete;

    /*!
     * Object name.
     */
    std::string d_object_name;

    /*!
     * Force magnitude.
     */
    double d_force_magnitude = 0.0;

    /*!
     * Force direction vector.
     */
    double d_force_direction[NDIM] = {0.0};

    /*!
     * Time-dependent parameters.
     */
    double d_start_time = 0.0;
    double d_ramp_time = 0.0;
};

#endif // #ifndef included_CustomForceFunction
