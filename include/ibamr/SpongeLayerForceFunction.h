// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
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

#ifndef included_IBAMR_SpongeLayerForceFunction
#define included_IBAMR_SpongeLayerForceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIArray.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAIDatabase.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAISideData.h>
#include <SAMRAIVariable.h>

#include <array>
#include <string>

namespace IBAMR
{
class INSHierarchyIntegrator;
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
 * \brief Class SpongeLayerForceFunction provides forcing at physical boundaries
 * that weakly imposes homogeneous Dirichlet boundary conditions.
 */
class SpongeLayerForceFunction : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    SpongeLayerForceFunction(const std::string& object_name,
                             SAMRAIPointer<SAMRAIDatabase> input_db,
                             const INSHierarchyIntegrator* fluid_solver,
                             SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry);

    /*!
     * \brief Destructor.
     */
    ~SpongeLayerForceFunction() = default;

    /*!
     * \name Methods to set the data.
     */
    //\{

    /*!
     * \note This concrete IBTK::CartGridFunction is time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAIPointer<SAMRAIVariable> var,
                        SAMRAIPointer<SAMRAIPatch> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAIPointer<SAMRAIPatchLevel> level = SAMRAIPointer<SAMRAIPatchLevel>(nullptr)) override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SpongeLayerForceFunction() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SpongeLayerForceFunction(const SpongeLayerForceFunction& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SpongeLayerForceFunction& operator=(const SpongeLayerForceFunction& that) = delete;

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatchCell(SAMRAIPointer<SAMRAICellData<double>> F_data,
                            SAMRAIPointer<SAMRAICellData<double>> U_current_data,
                            SAMRAIPointer<SAMRAICellData<double>> U_new_data,
                            double kappa,
                            SAMRAIPointer<SAMRAIPatch> patch);

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatchSide(SAMRAIPointer<SAMRAISideData<double>> F_data,
                            SAMRAIPointer<SAMRAISideData<double>> U_current_data,
                            SAMRAIPointer<SAMRAISideData<double>> U_new_data,
                            double kappa,
                            SAMRAIPointer<SAMRAIPatch> patch);

    std::array<SAMRAIArray<bool>, 2 * NDIM> d_forcing_enabled;
    std::array<double, 2 * NDIM> d_width;
    const INSHierarchyIntegrator* const d_fluid_solver;
    SAMRAIPointer<SAMRAICartesianGridGeometry> d_grid_geometry;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_SpongeLayerForceFunction
