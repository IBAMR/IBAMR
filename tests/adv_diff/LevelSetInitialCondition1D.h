// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_LevelSetInitialCondition1D
#define included_LevelSetInitialCondition1D

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibtk/muParserCartGridFunction.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class LevelSetInitialCondition1D provides an initial condition for
 * the level set function.
 */
class LevelSetInitialCondition1D : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    LevelSetInitialCondition1D(const std::string& object_name,
                               const IBTK::SAMRAIPointer<SAMRAI::geom::CartesianGridGeometryNd> grid_geom,
                               const IBTK::VectorNd& interface_loc,
                               const bool left_side = false);

    /*!
     * \brief Empty destructor.
     */
    ~LevelSetInitialCondition1D() = default;

    /*!
     * \brief Indicates whether the concrete LevelSetInitialCondition1D object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAIPointer<SAMRAI::hier::VariableNd> var,
                        SAMRAIPointer<PatchNd> patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAIPointer<PatchLevelNd> patch_level = SAMRAIPointer<PatchLevelNd>(NULL)) override;

    //\}

private:
    /*!
     * Deleted default constructor.
     */
    LevelSetInitialCondition1D() = delete;

    /*!
     * Deleted copy constructor.
     */
    LevelSetInitialCondition1D(const LevelSetInitialCondition1D& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    LevelSetInitialCondition1D& operator=(const LevelSetInitialCondition1D& that) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * The Cartesian grid geometry object provides the extents of the
     * computational domain.
     */
    IBTK::SAMRAIPointer<SAMRAI::geom::CartesianGridGeometryNd> d_grid_geom;

    /*!
     * Interface location.
     */
    IBTK::VectorNd d_interface_loc;

    /*!
     * Sepecify the side on which interface is located.
     */
    bool d_left_side = false;
};

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_LevelSetInitialCondition1D
