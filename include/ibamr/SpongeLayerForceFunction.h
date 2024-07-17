// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#include "ibtk/CartGridFunction.h"

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "PatchLevel.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"

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
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             const INSHierarchyIntegrator* fluid_solver,
                             SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometryNd> grid_geometry);

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
                        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableNd> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchNd> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevelNd> level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevelNd>(NULL)) override;

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
    void setDataOnPatchCell(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataNd<double> > F_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataNd<double> > U_current_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataNd<double> > U_new_data,
                            double kappa,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchNd> patch);

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatchSide(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideDataNd<double> > F_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideDataNd<double> > U_current_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideDataNd<double> > U_new_data,
                            double kappa,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchNd> patch);

    std::array<SAMRAI::tbox::Array<bool>, 2 * NDIM> d_forcing_enabled;
    std::array<double, 2 * NDIM> d_width;
    const INSHierarchyIntegrator* const d_fluid_solver;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometryNd> d_grid_geometry;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_SpongeLayerForceFunction
