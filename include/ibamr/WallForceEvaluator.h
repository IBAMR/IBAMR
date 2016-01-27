// Filename: WallForceEvaluator.h
// Created by Steven Delong
//
// Copyright (c) 2002-2014, Steven Delong, Boyce Griffith
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

#ifndef included_WallForceEvaluator
#define included_WallForceEvaluator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LNodeSetData.h"
#include "ibamr/IBLagrangianForceStrategy.h"
#include "ibamr/Wall.h"

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
                                IBTK::LDataManager* const l_data_manager);

private:
    // Default constructor, not implemented.
    WallForceEvaluator();

    // Copy constructor, not implemented.
    WallForceEvaluator(const WallForceEvaluator& from);

    // Assignment operator, not implemented.
    WallForceEvaluator& operator=(const WallForceEvaluator& that);

    // collection of walls:
    std::vector<Wall> d_walls_vec;

    // grid geometry, used when making walls:
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geometry;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_WallForceEvaluator
