// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2026 by the IBAMR developers
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

#ifndef included_IBAMR_VelocityCartGridFunction
#define included_IBAMR_VelocityCartGridFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/CartGridFunction.h>

#include <tbox/HDFDatabase.h>

#include <CartesianGridGeometry.h>

#include <string>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBAMR
{
/*!
 * \brief Class VelocityCartGridFunction provides an method to save a veloicity from an IBAMR simulation and use that
 * field for a transport solve.
 */
class VelocityCartGridFunction : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief The default constructor sets the name of the strategy object.
     */
    VelocityCartGridFunction(const std::string& object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM>> grid_geometry,
                             const int velocity_index,
                             const std::string& db_base_name);

    ~VelocityCartGridFunction() override = default;

    /*!
     * \brief Indicates whether the concrete CartGridFunction object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1) override;

    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(nullptr)) override;

    void saveVelocity(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                      int iteration_num,
                      double loop_time);

    double readVelocity(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy, int iteration_num);

    void setVelocityIndex(const int velocity_index);

private:
    std::string getIterationDir(int iteration_num) const;
    std::string getIterationFilename(int iteration_num) const;
    static void ensureDirectoryExists(const std::string& dir_path);
    static bool fileExists(const std::string& filename);

protected:
    std::string d_object_name;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM>> d_grid_geometry;
    int d_velocity_index;
    std::string d_db_base_name;
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_VelocityCartGridFunction
