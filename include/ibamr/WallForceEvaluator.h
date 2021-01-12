// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_WallForceEvaluator
#define included_IBAMR_WallForceEvaluator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBLagrangianForceStrategy.h"
#include "ibamr/Wall.h"

#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LNodeSetData.h"

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

#include <vector>

namespace IBTK
{
class LData;
class LDataManager;
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

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
class WallForceEvaluator : public IBLagrangianForceStrategy
{
public:
    // Constructor.
    WallForceEvaluator(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                       SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry);

    // Register a WallForceFcnPtr to be used by all walls.
    void registerWallForceFcn(Wall::WallForceFcnPtr wall_force_fcn);

    // add a wall to the WallForceEvaluator
    void addWall(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> wall_db, double wall_ghost_dist);

    // compute forces from all walls
    void computeLagrangianForce(SAMRAI::tbox::Pointer<IBTK::LData> F_data,
                                SAMRAI::tbox::Pointer<IBTK::LData> X_data,
                                SAMRAI::tbox::Pointer<IBTK::LData> U_data,
                                const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                const int level_number,
                                const double data_time,
                                IBTK::LDataManager* const l_data_manager) override;

private:
    // Default constructor, not implemented.
    WallForceEvaluator() = delete;

    // Copy constructor, not implemented.
    WallForceEvaluator(const WallForceEvaluator& from) = delete;

    // Assignment operator, not implemented.
    WallForceEvaluator& operator=(const WallForceEvaluator& that) = delete;

    // collection of walls:
    std::vector<Wall> d_walls_vec;

    // grid geometry, used when making walls:
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geometry;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_WallForceEvaluator
