// Filename: SimplifiedIBFEMethod.C
// Created on 11 Sep 2012 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
// POSSIBILITY OF SUCH DAMAGE.

#include "SimplifiedIBFEMethod.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LEInteractor.h>

// LIBMESH INCLUDES
#include <boundary_info.h>
#include <dense_vector.h>
#include <fe_interface.h>
#include <mesh.h>
#include <petsc_vector.h>
#include <string_to_enum.h>
using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of SimplifiedIBFEMethod restart file data.
static const int SIMPLIFIED_IBFE_METHOD_VERSION = 1;

inline
short int
get_dirichlet_bdry_ids(
    const std::vector<short int>& bdry_ids)
{
    short int dirichlet_bdry_ids = 0;
    for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
    {
        const short int bdry_id = *cit;
        if      (bdry_id == FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID  ) dirichlet_bdry_ids |= FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID;
        else if (bdry_id == FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID  ) dirichlet_bdry_ids |= FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID;
        else if (bdry_id == FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID  ) dirichlet_bdry_ids |= FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID;
        else if (bdry_id == FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID ) dirichlet_bdry_ids |= FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID;
        else if (bdry_id == FEDataManager::ZERO_DISPLACEMENT_XZ_BDRY_ID ) dirichlet_bdry_ids |= FEDataManager::ZERO_DISPLACEMENT_XZ_BDRY_ID;
        else if (bdry_id == FEDataManager::ZERO_DISPLACEMENT_YZ_BDRY_ID ) dirichlet_bdry_ids |= FEDataManager::ZERO_DISPLACEMENT_YZ_BDRY_ID;
        else if (bdry_id == FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID) dirichlet_bdry_ids |= FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID;
    }
    return dirichlet_bdry_ids;
}// get_dirichlet_bdry_ids

static const double MIN_POINTS = 2.0;
static const double POINT_FACTOR = 2.0;

inline double
get_elem_hmax(
    Elem* elem,
    const blitz::Array<double,2>& X_node)
{
    static const int MAX_NODES = (NDIM == 2 ? 9 : 27);
    Point s_node_cache[MAX_NODES];
    const int n_node = elem->n_nodes();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(n_node <= MAX_NODES);
#endif
    for (int k = 0; k < n_node; ++k)
    {
        s_node_cache[k] = elem->point(k);
        Point& X = elem->point(k);
        for (int d = 0; d < NDIM; ++d)
        {
            X(d) = X_node(k,d);
        }
    }
    const double hmax = elem->hmax();
    for (int k = 0; k < n_node; ++k)
    {
        elem->point(k) = s_node_cache[k];
    }
    return hmax;
}// get_elem_hmax

#if 0
static const int kernel_width = 2;

inline double
kernel(
    double x)
{
    x += 2.;
    const double x2 = x*x;
    const double x3 = x*x2;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .1666666666666667*x3;
    else if (x <= 2.)
        return 2.*x2-.5000000000000000*x3-2.*x+.6666666666666667;
    else if (x <= 3.)
        return 10.*x-4.*x2+.5000000000000000*x3-7.333333333333333;
    else if (x <= 4.)
        return 10.66666666666667-8.*x+2.*x2-.1666666666666667*x3;
    else
        return 0.;
}// kernel

inline double
kernel_diff(
    double x)
{
    x += 2.;
    const double x2 = x*x;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .5000000000000000*x2;
    else if (x <= 2.)
        return 4.*x-1.500000000000000*x2-2.;
    else if (x <= 3.)
        return 10.-8.*x+1.500000000000000*x2;
    else if (x <= 4.)
        return -8.+4.*x-.5000000000000000*x2;
    else
        return 0.;
}// kernel_diff
#endif

#if 0
static const int kernel_width = 3;

inline double
kernel(
    double x)
{
    x += 3.;
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    const double x5 = x*x4;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .8333333333333333e-2*x5;
    else if (x <= 2.)
        return .2500000000000000*x4-.4166666666666667e-1*x5-.5000000000000000*x3+.5000000000000000*x2-.2500000000000000*x+.5000000000000000e-1;
    else if (x <= 3.)
        return 4.500000000000000*x3-1.*x4+.8333333333333333e-1*x5-9.500000000000000*x2+9.750000000000000*x-3.950000000000000;
    else if (x <= 4.)
        return 35.50000000000000*x2-10.50000000000000*x3+1.500000000000000*x4-.8333333333333333e-1*x5-57.75000000000000*x+36.55000000000000;
    else if (x <= 5.)
        return 102.2500000000000*x-44.50000000000000*x2+9.500000000000000*x3-1.*x4+.4166666666666667e-1*x5-91.45000000000000;
    else if (x <= 6.)
        return 64.80000000000000-54.*x+18.*x2-3.*x3+.2500000000000000*x4-.8333333333333333e-2*x5;
    else
        return 0.;
}// kernel

inline double
kernel_diff(
    double x)
{
    x += 3.;
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    if (x < 0.)
        return 0.;
    else if (x <= 1.)
        return .4166666666666667e-1*x4;
    else if (x <= 2.)
        return -.2500000000000000+x-1.500000000000000*x2+x3-.2083333333333333*x4;
    else if (x <= 3.)
        return 9.750000000000000-19.*x+13.50000000000000*x2-4.*x3+.4166666666666667*x4;
    else if (x <= 4.)
        return -57.75000000000000+71.*x-31.50000000000000*x2+6.*x3-.4166666666666667*x4;
    else if (x <= 5.)
        return 102.2500000000000-89.*x+28.50000000000000*x2-4.*x3+.2083333333333333*x4;
    else if (x <= 6.)
        return -54.+36.*x-9.*x2+x3-.4166666666666667e-1*x4;
    else
        return 0.;
}// kernel_diff
#endif

#if 1
static const int kernel_width = 4;

inline double
kernel(
    double x)
{
    x += 4.;
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    const double x5 = x*x4;
    const double x6 = x*x5;
    const double x7 = x*x6;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .1984126984126984e-3*x7;
    else if (x <= 2.)
        return .1111111111111111e-1*x6-.1388888888888889e-2*x7-.3333333333333333e-1*x5+.5555555555555556e-1*x4-.5555555555555556e-1*x3+.3333333333333333e-1*x2-.1111111111111111e-1*x+.1587301587301587e-2;
    else if (x <= 3.)
        return .4333333333333333*x5-.6666666666666667e-1*x6+.4166666666666667e-2*x7-1.500000000000000*x4+3.055555555555556*x3-3.700000000000000*x2+2.477777777777778*x-.7095238095238095;
    else if (x <= 4.)
        return 9.*x4-1.666666666666667*x5+.1666666666666667*x6-.6944444444444444e-2*x7-28.44444444444444*x3+53.*x2-54.22222222222222*x+23.59047619047619;
    else if (x <= 5.)
        return 96.*x3-22.11111111111111*x4+3.*x5-.2222222222222222*x6+.6944444444444444e-2*x7-245.6666666666667*x2+344.*x-203.9650793650794;
    else if (x <= 6.)
        return 483.5000000000000*x2-147.0555555555556*x3+26.50000000000000*x4-2.833333333333333*x5+.1666666666666667*x6-.4166666666666667e-2*x7-871.2777777777778*x+664.0904761904762;
    else if (x <= 7.)
        return 943.1222222222222*x-423.7000000000000*x2+104.9444444444444*x3-15.50000000000000*x4+1.366666666666667*x5-.6666666666666667e-1*x6+.1388888888888889e-2*x7-891.1095238095238;
    else if (x <= 8.)
        return 416.1015873015873-364.0888888888889*x+136.5333333333333*x2-28.44444444444444*x3+3.555555555555556*x4-.2666666666666667*x5+.1111111111111111e-1*x6-.1984126984126984e-3*x7;
    else
        return 0.;
}// kernel

inline double
kernel_diff(
    double x)
{
    x += 4.;
    const double x2 = x*x;
    const double x3 = x*x2;
    const double x4 = x*x3;
    const double x5 = x*x4;
    const double x6 = x*x5;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .1388888888888889e-2*x6;
    else if (x <= 2.)
        return -.1111111111111111e-1+.6666666666666667e-1*x-.1666666666666667*x2+.2222222222222222*x3-.1666666666666667*x4+.6666666666666667e-1*x5-.9722222222222222e-2*x6;
    else if (x <= 3.)
        return 2.477777777777778-7.400000000000000*x+9.166666666666667*x2-6.*x3-.4000000000000000*x5+.2916666666666667e-1*x6+2.166666666666667*x4;
    else if (x <= 4.)
        return -54.22222222222222+106.*x-85.33333333333333*x2+36.*x3-8.333333333333333*x4+x5-.4861111111111111e-1*x6;
    else if (x <= 5.)
        return 344.-491.3333333333333*x+288.*x2-88.44444444444444*x3+15.*x4-1.333333333333333*x5+.4861111111111111e-1*x6;
    else if (x <= 6.)
        return -871.2777777777778+967.*x-441.1666666666667*x2+106.*x3-14.16666666666667*x4+x5-.2916666666666667e-1*x6;
    else if (x <= 7.)
        return 943.1222222222222-847.4000000000000*x+314.8333333333333*x2-62.*x3+6.833333333333333*x4-.4000000000000000*x5+.9722222222222222e-2*x6;
    else if (x <= 8.)
        return -364.0888888888889+273.0666666666667*x-85.33333333333333*x2+14.22222222222222*x3-1.333333333333333*x4+.6666666666666667e-1*x5-.1388888888888889e-2*x6;
    else
        return 0.;
}// kernel_diff
#endif

}

const std::string SimplifiedIBFEMethod::       COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string SimplifiedIBFEMethod::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string SimplifiedIBFEMethod::     VELOCITY_SYSTEM_NAME = "IB velocity system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

SimplifiedIBFEMethod::SimplifiedIBFEMethod(
    const std::string& object_name,
    Pointer<Database> input_db,
    Mesh* mesh,
    int max_level_number,
    bool register_for_restart)
    : d_num_parts(1)
{
    commonConstructor(object_name, input_db, std::vector<Mesh*>(1,mesh), max_level_number, register_for_restart);
    return;
}// SimplifiedIBFEMethod

SimplifiedIBFEMethod::SimplifiedIBFEMethod(
    const std::string& object_name,
    Pointer<Database> input_db,
    const std::vector<Mesh*>& meshes,
    int max_level_number,
    bool register_for_restart)
    : d_num_parts(meshes.size())
{
    commonConstructor(object_name, input_db, meshes, max_level_number, register_for_restart);
    return;
}// SimplifiedIBFEMethod

SimplifiedIBFEMethod::~SimplifiedIBFEMethod()
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        delete d_equation_systems[part];
    }
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
}// ~SimplifiedIBFEMethod

FEDataManager*
SimplifiedIBFEMethod::getFEDataManager(
    const unsigned int part) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(part < d_num_parts);
#endif
    return d_fe_data_managers[part];
}// getFEDataManager

void
SimplifiedIBFEMethod::registerInitialCoordinateMappingFunction(
    CoordinateMappingFcnPtr coordinate_mapping_fcn,
    void* coordinate_mapping_fcn_ctx,
    const unsigned int part)
{
    d_coordinate_mapping_fcns    [part] = coordinate_mapping_fcn;
    d_coordinate_mapping_fcn_ctxs[part] = coordinate_mapping_fcn_ctx;
    return;
}// registerInitialCoordinateMappingFunction

void
SimplifiedIBFEMethod::registerPK1StressTensorFunction(
    PK1StressFcnPtr PK1_stress_fcn,
    std::vector<unsigned int> PK1_stress_fcn_systems,
    void* PK1_stress_fcn_ctx,
    const unsigned int part)
{
    d_PK1_stress_fcns       [part] = PK1_stress_fcn;
    d_PK1_stress_fcn_systems[part] = PK1_stress_fcn_systems;
    d_PK1_stress_fcn_ctxs   [part] = PK1_stress_fcn_ctx;
    return;
}// registerPK1StressTensorFunction

const IntVector<NDIM>&
SimplifiedIBFEMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
}// getMinimumGhostCellWidth

void
SimplifiedIBFEMethod::setupTagBuffer(
    Array<int>& tag_buffer,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    const int finest_hier_ln = gridding_alg->getMaxLevels()-1;
    const int tsize = tag_buffer.size();
    tag_buffer.resizeArray(finest_hier_ln);
    for (int i = tsize; i < finest_hier_ln; ++i) tag_buffer[i] = 0;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const int gcw    = d_fe_data_managers[part]->getGhostCellWidth().max();
        const int tag_ln = d_fe_data_managers[part]->getLevelNumber()-1;
        if (tag_ln >= 0 && tag_ln < finest_hier_ln)
        {
            tag_buffer[tag_ln] = std::max(tag_buffer[tag_ln], gcw);
        }
    }
    for (int ln = finest_hier_ln-2; ln >= 0; --ln)
    {
        tag_buffer[ln] = std::max(tag_buffer[ln], tag_buffer[ln+1]/gridding_alg->getRatioToCoarserLevel(ln+1).max()+1);
    }
    return;
}// setupTagBuffer

void
SimplifiedIBFEMethod::preprocessIntegrateData(
    double current_time,
    double new_time,
    int /*num_cycles*/)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time+0.5*(new_time-current_time);

    // Extract the FE data.
    d_X_systems      .resize(d_num_parts);
    d_X_current_vecs .resize(d_num_parts);
    d_X_new_vecs     .resize(d_num_parts);
    d_X_half_vecs    .resize(d_num_parts);
    d_X_IB_ghost_vecs.resize(d_num_parts);
    d_U_systems      .resize(d_num_parts);
    d_U_current_vecs .resize(d_num_parts);
    d_U_new_vecs     .resize(d_num_parts);
    d_U_half_vecs    .resize(d_num_parts);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_X_systems      [part] = &d_equation_systems[part]->get_system(  COORDS_SYSTEM_NAME);
        d_X_current_vecs [part] = dynamic_cast<PetscVector<double>*>(d_X_systems       [part]->solution.get());
        d_X_new_vecs     [part] = dynamic_cast<PetscVector<double>*>(d_X_current_vecs  [part]->clone().release());  // WARNING: must be manually deleted
        d_X_half_vecs    [part] = dynamic_cast<PetscVector<double>*>(d_X_systems       [part]->current_local_solution.get());
        d_X_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(d_fe_data_managers[part]->buildGhostedCoordsVector());
        d_U_systems      [part] = &d_equation_systems[part]->get_system(VELOCITY_SYSTEM_NAME);
        d_U_current_vecs [part] = dynamic_cast<PetscVector<double>*>(d_U_systems       [part]->solution.get());
        d_U_new_vecs     [part] = dynamic_cast<PetscVector<double>*>(d_U_current_vecs  [part]->clone().release());  // WARNING: must be manually deleted
        d_U_half_vecs    [part] = dynamic_cast<PetscVector<double>*>(d_U_systems       [part]->current_local_solution.get());

        // Initialize X^{n+1/2} and X^{n+1} to equal X^{n}, and initialize
        // U^{n+1/2} and U^{n+1} to equal U^{n}.
        d_X_current_vecs[part]->localize(*d_X_half_vecs[part]);
        d_X_half_vecs[part]->close();
        d_X_current_vecs[part]->localize(*d_X_new_vecs[part]);
        d_X_new_vecs[part]->close();
        d_U_current_vecs[part]->localize(*d_U_half_vecs[part]);
        d_U_half_vecs[part]->close();
        d_U_current_vecs[part]->localize(*d_U_new_vecs[part]);
        d_U_new_vecs[part]->close();
    }
    return;
}// preprocessIntegrateData

void
SimplifiedIBFEMethod::postprocessIntegrateData(
    double /*current_time*/,
    double /*new_time*/,
    int /*num_cycles*/)
{
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        // Reset time-dependent Lagrangian data.
        (*d_X_current_vecs[part]) = (*d_X_new_vecs[part]);
        (*d_U_current_vecs[part]) = (*d_U_new_vecs[part]);

        d_X_systems[part]->solution->localize(*d_X_systems[part]->current_local_solution);
        d_U_systems[part]->solution->localize(*d_U_systems[part]->current_local_solution);

        // Update the coordinate mapping dX = X - s.
        updateCoordinateMapping(part);

        // Deallocate Lagrangian scratch data.
        delete d_X_new_vecs[part];
        delete d_U_new_vecs[part];
    }
    d_X_systems      .clear();
    d_X_current_vecs .clear();
    d_X_new_vecs     .clear();
    d_X_half_vecs    .clear();
    d_X_IB_ghost_vecs.clear();
    d_U_systems      .clear();
    d_U_current_vecs .clear();
    d_U_new_vecs     .clear();
    d_U_half_vecs    .clear();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time     = std::numeric_limits<double>::quiet_NaN();
    d_half_time    = std::numeric_limits<double>::quiet_NaN();
    return;
}// postprocessIntegrateData

void
SimplifiedIBFEMethod::interpolateVelocity(
    const int u_data_idx,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
    const double data_time)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* X_vec = NULL;
        PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        PetscVector<double>* U_vec = NULL;
        if (MathUtilities<double>::equalEps(data_time, d_current_time))
        {
            X_vec = d_X_current_vecs[part];
            U_vec = d_U_current_vecs[part];
        }
        else if (MathUtilities<double>::equalEps(data_time, d_half_time))
        {
            X_vec = d_X_half_vecs[part];
            U_vec = d_U_half_vecs[part];
        }
        else if (MathUtilities<double>::equalEps(data_time, d_new_time))
        {
            X_vec = d_X_new_vecs[part];
            U_vec = d_U_new_vecs[part];
        }
        X_vec->localize(*X_ghost_vec);
        System& X_system = d_equation_systems[part]->get_system(COORDS_SYSTEM_NAME);
        X_system.get_dof_map().enforce_constraints_exactly(X_system, X_ghost_vec);
        X_ghost_vec->close();
        projectMaterialVelocity(u_data_idx, *U_vec, *X_ghost_vec, data_time, part);
    }
    return;
}// interpolateVelocity

void
SimplifiedIBFEMethod::eulerStep(
    const double current_time,
    const double new_time)
{
    const double dt = new_time-current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec()); IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());  IBTK_CHKERRQ(ierr);
        d_X_new_vecs [part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
}// eulerStep

void
SimplifiedIBFEMethod::midpointStep(
    const double current_time,
    const double new_time)
{
    const double dt = new_time-current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_half_vecs[part]->vec(), d_X_current_vecs[part]->vec()); IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec()); IBTK_CHKERRQ(ierr);
        d_X_new_vecs [part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
}// midpointStep

void
SimplifiedIBFEMethod::trapezoidalStep(
    const double current_time,
    const double new_time)
{
    const double dt = new_time-current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), 0.5*dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec()); IBTK_CHKERRQ(ierr);
        ierr = VecAXPY( d_X_new_vecs[part]->vec(), 0.5*dt, d_U_new_vecs    [part]->vec()); IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec()); IBTK_CHKERRQ(ierr);
        d_X_new_vecs [part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
}// trapezoidalStep

void
SimplifiedIBFEMethod::computeLagrangianForce(
    const double /*data_time*/)
{
    // intentionally blank
    return;
}// computeLagrangianForce

void
SimplifiedIBFEMethod::spreadForce(
    const int f_data_idx,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
    const double data_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
#endif
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Get the IB ghosted position vector.
        PetscVector<double>& X_ghost_vec = *d_X_IB_ghost_vecs[part];
        d_X_half_vecs[part]->localize(X_ghost_vec);
        System& X_system = d_equation_systems[part]->get_system(COORDS_SYSTEM_NAME);
        X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_ghost_vec);
        X_ghost_vec.close();
        projectTotalForceDensity(f_data_idx, X_ghost_vec, data_time, part);
    }
    return;
}// spreadForce

void
SimplifiedIBFEMethod::initializeFEData()
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Initialize FE equation systems.
        EquationSystems* equation_systems = d_equation_systems[part];
        equation_systems->init();
        initializeCoordinates(part);
        updateCoordinateMapping(part);

        // Assemble systems.
        System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
        X_system.assemble_before_solve = false;
        X_system.assemble();

        System& U_system = equation_systems->get_system<System>(VELOCITY_SYSTEM_NAME);
        U_system.assemble_before_solve = false;
        U_system.assemble();

        System& X_mapping_system = equation_systems->get_system<System>(COORD_MAPPING_SYSTEM_NAME);
        X_mapping_system.assemble_before_solve = false;
        X_mapping_system.assemble();

        // Set up boundary conditions.  Specifically, add appropriate boundary
        // IDs to the BoundaryInfo object associated with the mesh, and add DOF
        // constraints for the nodal forces and velocities.
        const MeshBase& mesh = equation_systems->get_mesh();
        DofMap& U_dof_map = U_system.get_dof_map();
        const unsigned int U_sys_num = U_system.number();
        MeshBase::const_element_iterator       el_it  = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for ( ; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = elem->neighbor(side) == NULL;
                if (!at_mesh_bdry) continue;

                static const short int dirichlet_bdry_id_set[3] = { FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID , FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID , FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID };
                const short int dirichlet_bdry_ids = get_dirichlet_bdry_ids(mesh.boundary_info->boundary_ids(elem, side));
                if (!dirichlet_bdry_ids) continue;

                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (!elem->is_node_on_side(n, side)) continue;

                    Node* node = elem->get_node(n);
                    mesh.boundary_info->add_node(node, dirichlet_bdry_ids);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (!(dirichlet_bdry_ids & dirichlet_bdry_id_set[d])) continue;
                        if (node->n_dofs(U_sys_num) > 0)
                        {
                            const int U_dof_index = node->dof_number(U_sys_num,d,0);
                            DofConstraintRow U_constraint_row;
                            U_constraint_row[U_dof_index] = 1.0;
                            U_dof_map.add_constraint_row(U_dof_index, U_constraint_row, 0.0, false);
                        }
                    }
                }
            }
        }
    }

    d_is_initialized = true;
    return;
}// initializeFEData

void
SimplifiedIBFEMethod::initializePatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
    int /*u_data_idx*/,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
    int /*integrator_step*/,
    double /*init_data_time*/,
    bool /*initial_time*/)
{
    // Cache pointers to the patch hierarchy and gridding algorithm.
    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Initialize the FE data manager.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->reinitElementMappings();
    }

    d_is_initialized = true;
    return;
}// initializePatchHierarchy

void
SimplifiedIBFEMethod::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > load_balancer,
    int workload_data_idx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(load_balancer);
#endif
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->registerLoadBalancer(load_balancer, workload_data_idx);
    }
    return;
}// registerLoadBalancer

void
SimplifiedIBFEMethod::updateWorkloadEstimates(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    int /*workload_data_idx*/)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->updateWorkloadEstimates();
    }
    return;
}// updateWorkloadEstimates

void
SimplifiedIBFEMethod::beginDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
}// beginDataRedistribution

void
SimplifiedIBFEMethod::endDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    if (d_is_initialized)
    {
        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            d_fe_data_managers[part]->reinitElementMappings();
        }
    }
    return;
}// endDataRedistribution

void
SimplifiedIBFEMethod::initializeLevelData(
    Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    int level_number,
    double init_data_time,
    bool can_be_refined,
    bool initial_time,
    Pointer<BasePatchLevel<NDIM> > old_level,
    bool allocate_data)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
        d_fe_data_managers[part]->resetLevels(0,finest_hier_level);
        d_fe_data_managers[part]->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
        if (d_load_balancer && level_number == d_fe_data_managers[part]->getLevelNumber())
        {
            d_load_balancer->setWorkloadPatchDataIndex(d_workload_idx, level_number);
            d_fe_data_managers[part]->updateWorkloadEstimates(level_number, level_number);
        }
    }
    return;
}// initializeLevelData

void
SimplifiedIBFEMethod::resetHierarchyConfiguration(
    Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    int coarsest_level,
    int /*finest_level*/)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
        d_fe_data_managers[part]->resetLevels(0,hierarchy->getFinestLevelNumber());
        d_fe_data_managers[part]->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_hier_level);
    }
    return;
}// resetHierarchyConfiguration

void
SimplifiedIBFEMethod::applyGradientDetector(
    Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    int level_number,
    double error_data_time,
    int tag_index,
    bool initial_time,
    bool uses_richardson_extrapolation_too)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }
    return;
}// applyGradientDetector

void
SimplifiedIBFEMethod::putToDatabase(
    Pointer<Database> db)
{
    db->putInteger("SIMPLIFIED_IBFE_METHOD_VERSION", SIMPLIFIED_IBFE_METHOD_VERSION);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    db->putString("d_fe_family", Utility::enum_to_string<FEFamily>(d_fe_family));
    db->putString("d_fe_order", Utility::enum_to_string<Order>(d_fe_order));
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

void
SimplifiedIBFEMethod::projectMaterialVelocity(
    const int u_data_idx,
    PetscVector<double>& U_vec,
    PetscVector<double>& X_ghost_vec,
    const double /*time*/,
    const unsigned int part)
{
    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const int dim = mesh.mesh_dimension();
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, FIRST);

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& system = equation_systems->get_system(VELOCITY_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(dof_map.variable_type(d) == dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) dof_indices(d).reserve(27);
    AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();

#ifdef DEBUG_CHECK_ASSERTIONS
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(dof_map.variable_type(d) == X_dof_map.variable_type(d));
#endif

    // Loop over the patches to comptue the righ-hand side term
    AutoPtr<NumericVector<double> > U_rhs_vec = U_vec.zero_clone();
    std::vector<DenseVector<double> > U_rhs_e(NDIM);
    const blitz::Array<blitz::Array<Elem*,1>,1>& active_patch_element_map = d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    VectorValue<double> X_cell, X_qp;
    blitz::Array<double,2> X_node;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const blitz::Array<Elem*,1>& patch_elems = active_patch_element_map(local_patch_num);
        const unsigned int num_active_patch_elems = patch_elems.size();
        if (!num_active_patch_elems) continue;

        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<SideData<NDIM,double> > u_data = patch->getPatchData(u_data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const x_lower = patch_geom->getXLower();
        const double* const x_upper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();
        const double dx_min = *std::min_element(dx,dx+NDIM);
        Box<NDIM> side_boxes[NDIM];
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            side_boxes[axis] = patch_box;
            if (patch_geom->getTouchesRegularBoundary(axis, /*upper*/ 1)) side_boxes[axis].growUpper(axis,1);
        }
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems(e_idx);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dof_map.dof_indices(elem, dof_indices(d), d);
                U_rhs_e[d].resize(dof_indices(d).size());  // NOTE: DenseVector::resize() automatically zeroes the vector contents.
            }
            get_values_for_interpolation(X_node, X_ghost_vec, dof_indices);
            const double hmax = get_elem_hmax(elem, X_node);
            const int npts = std::max(MIN_POINTS, std::ceil(POINT_FACTOR*hmax/dx_min));
            const Order order = static_cast<Order>(std::min(2*npts-1,static_cast<int>(FORTYTHIRD)));
            if (order != qrule->get_order())
            {
                qrule = QBase::build(QGAUSS, dim, order);
                fe->attach_quadrature_rule(qrule.get());
            }
            fe->reinit(elem);
            for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
            {
                interpolate(X_qp,qp,X_node,phi);

                // Weight u using a smooth kernel function evaluated about X_qp.
                //
                // WARNING: As written here, this implicitly imposes u = 0 in
                // the ghost cell region.
                const Index<NDIM> i = IndexUtilities::getCellIndex(&X_qp(0), x_lower, x_upper, dx, patch_box.lower(), patch_box.upper());
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X_cell(d) = x_lower[d] + dx[d]*(static_cast<double>(i(d)-patch_box.lower(d))+0.5);
                }
                VectorValue<double> U;
                for (unsigned int component = 0; component < NDIM; ++component)
                {
                    blitz::Array<double,1> phi[NDIM];
                    Box<NDIM> box(i,i);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == component)
                        {
                            box.lower(d) -= (kernel_width-1);
                            box.upper(d) += (kernel_width  );
                        }
                        else
                        {
                            box.lower(d) -= (X_qp(d) <= X_cell(d) ? kernel_width : kernel_width-1);
                            box.upper(d) += (X_qp(d) >= X_cell(d) ? kernel_width : kernel_width-1);
                        }
                        phi[d].resize(blitz::Range(box.lower(d),box.upper(d)));
                    }
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        for (int i = box.lower(d); i <= box.upper(d); ++i)
                        {
                            const double x_grid = x_lower[d] + dx[d]*(static_cast<double>(i-patch_box.lower(d))+(d == component ? 0.0 : 0.5));
                            const double del = x_grid - X_qp(d);
                            phi[d](i) = kernel(del/dx[d]);
                        }
                    }
                    for (Box<NDIM>::Iterator b(box*side_boxes[component]); b; b++)
                    {
                        const Index<NDIM>& i = b();
                        double w = 1.0;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            w *= phi[d](i(d));
                        }
                        U(component) += (*u_data)(SideIndex<NDIM>(i, component, SideIndex<NDIM>::Lower))*w;
                    }
                }
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    for (unsigned int k = 0; k < dof_indices(d).size(); ++k)
                    {
                        U_rhs_e[d](k) += U(d)*phi[k][qp]*JxW[qp];
                    }
                }
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dof_map.constrain_element_vector(U_rhs_e[d], dof_indices(d));
                U_rhs_vec->add_vector(U_rhs_e[d], dof_indices(d));
            }
        }
    }

    // Solve for the nodal values.
    d_fe_data_managers[part]->computeL2Projection(U_vec, *U_rhs_vec, COORDS_SYSTEM_NAME, d_use_consistent_mass_matrix);
    return;
}// projectMaterialVelocity

void
SimplifiedIBFEMethod::projectTotalForceDensity(
    const int f_data_idx,
    PetscVector<double>& X_ghost_vec,
    const double time,
    const unsigned int part)
{
    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const int dim = mesh.mesh_dimension();
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, FIRST);

    // Setup extra data needed to compute stresses/forces.
    std::vector<NumericVector<double>*> PK1_stress_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_PK1_stress_fcn_systems[part].begin(); cit != d_PK1_stress_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        PK1_stress_fcn_data.push_back(d_fe_data_managers[part]->buildGhostedSolutionVector(system.name()));
    }

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d) TBOX_ASSERT(dof_map.variable_type(d) == dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> dof_indices(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) dof_indices(d).reserve(27);
    AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<Point>& q_point = fe->get_xyz();
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi = fe->get_dphi();

    // Loop over the patches to comptue the forces on the grid.
    const blitz::Array<blitz::Array<Elem*,1>,1>& active_patch_element_map = d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    TensorValue<double> PP, FF, tau;
    VectorValue<double> X_cell, X_qp;
    blitz::Array<double,2> X_node;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const blitz::Array<Elem*,1>& patch_elems = active_patch_element_map(local_patch_num);
        const unsigned int num_active_patch_elems = patch_elems.size();
        if (!num_active_patch_elems) continue;

        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        Box<NDIM> side_boxes[NDIM];
        for (unsigned int d = 0; d < NDIM; ++d) side_boxes[d] = SideGeometry<NDIM>::toSideBox(patch_box,d);
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const x_lower = patch_geom->getXLower();
        const double* const x_upper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();
        const double dx_min = *std::min_element(dx,dx+NDIM);
        double dV_c = 1.0; for (unsigned int d = 0; d < NDIM; ++d) dV_c *= dx[d];
        for (unsigned int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems(e_idx);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dof_map.dof_indices(elem, dof_indices(d), d);
            }
            get_values_for_interpolation(X_node, X_ghost_vec, dof_indices);
            const double hmax = get_elem_hmax(elem, X_node);
            const int npts = std::max(MIN_POINTS, std::ceil(POINT_FACTOR*hmax/dx_min));
            const Order order = static_cast<Order>(std::min(2*npts-1,static_cast<int>(FORTYTHIRD)));
            if (order != qrule->get_order())
            {
                qrule = QBase::build(QGAUSS, dim, order);
                fe->attach_quadrature_rule(qrule.get());
            }
            fe->reinit(elem);
            for (unsigned int qp = 0; qp < qrule->n_points(); ++qp)
            {
                interpolate(X_qp,qp,X_node,phi);
                const Point& s_qp = q_point[qp];
                interpolate(X_qp,qp,X_node,phi);
                jacobian(FF,qp,X_node,dphi);
                if (d_PK1_stress_fcns[part])
                {
                    d_PK1_stress_fcns[part](PP,FF,X_qp,s_qp,elem,X_ghost_vec,PK1_stress_fcn_data,time,d_PK1_stress_fcn_ctxs[part]);
                }
                else PP.zero();
                tau = PP*FF.transpose();

                // Weight tau using a smooth kernel function evaluated about
                // X_qp.
                const Index<NDIM> i = IndexUtilities::getCellIndex(&X_qp(0), x_lower, x_upper, dx, patch_box.lower(), patch_box.upper());
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X_cell(d) = x_lower[d] + dx[d]*(static_cast<double>(i(d)-patch_box.lower(d))+0.5);
                }
                for (unsigned int component = 0; component < NDIM; ++component)
                {
                    blitz::Array<double,1> phi[NDIM], dphi[NDIM];
                    Box<NDIM> box(i,i);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == component)
                        {
                            box.lower(d) -= (kernel_width-1);
                            box.upper(d) += (kernel_width  );
                        }
                        else
                        {
                            box.lower(d) -= (X_qp(d) <= X_cell(d) ? kernel_width : kernel_width-1);
                            box.upper(d) += (X_qp(d) >= X_cell(d) ? kernel_width : kernel_width-1);
                        }
                        phi [d].resize(blitz::Range(box.lower(d),box.upper(d)));
                        dphi[d].resize(blitz::Range(box.lower(d),box.upper(d)));
                    }
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        for (int i = box.lower(d); i <= box.upper(d); ++i)
                        {
                            const double x_grid = x_lower[d] + dx[d]*(static_cast<double>(i-patch_box.lower(d))+(d == component ? 0.0 : 0.5));
                            const double del = x_grid - X_qp(d);
                            phi [d](i) =  kernel(del/dx[d]);
                            dphi[d](i) = -kernel_diff(del/dx[d])/dx[d];  // note the minus sign is because we are differentiating wrt X
                        }
                    }
                    for (Box<NDIM>::Iterator b(box*side_boxes[component]); b; b++)
                    {
                        const Index<NDIM>& i = b();
                        double f = 0.0;
                        for (unsigned int k = 0; k < NDIM; ++k)
                        {
                            double dw_dx_k = 1.0;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                if (d == k) dw_dx_k *= dphi[d](i(d));
                                else        dw_dx_k *=  phi[d](i(d));
                            }
                            f -= tau(component,k)*dw_dx_k;
                        }
                        (*f_data)(SideIndex<NDIM>(i, component, SideIndex<NDIM>::Lower)) += f*JxW[qp]/dV_c;
                    }
                }
            }
        }
    }
    return;
}// projectTotalForceDensity

void
SimplifiedIBFEMethod::initializeCoordinates(
    const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num) > 0)
        {
            libmesh_assert(n->n_vars(X_sys_num) == NDIM);
            const Point& s = *n;
            Point X = s;
            if (d_coordinate_mapping_fcns[part])
            {
                d_coordinate_mapping_fcns[part](X, s, d_coordinate_mapping_fcn_ctxs[part]);
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(X_sys_num,d,0);
                X_coords.set(dof_index,X(d));
            }
        }
    }
    X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
    X_coords.close();
    return;
}// initializeCoordinates

void
SimplifiedIBFEMethod::updateCoordinateMapping(
    const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    System& X_mapping_system = equation_systems->get_system(COORD_MAPPING_SYSTEM_NAME);
    const unsigned int X_mapping_sys_num = X_mapping_system.number();
    NumericVector<double>& dX_coords = *X_mapping_system.solution;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num) > 0)
        {
            libmesh_assert(n->n_vars(X_sys_num) == NDIM);
            libmesh_assert(n->n_vars(X_mapping_sys_num) == NDIM);
            const Point& s = *n;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int X_dof_index = n->dof_number(X_sys_num,d,0);
                const int dX_dof_index = n->dof_number(X_mapping_sys_num,d,0);
                dX_coords.set(dX_dof_index,X_coords(X_dof_index)-s(d));
            }
        }
    }
    dX_coords.close();
    return;
}// updateCoordinateMapping

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SimplifiedIBFEMethod::commonConstructor(
    const std::string& object_name,
    Pointer<Database> input_db,
    const std::vector<libMesh::Mesh*>& meshes,
    int max_level_number,
    bool register_for_restart)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }

    // Set some default values.
    d_ghosts = 0;
    d_use_consistent_mass_matrix = true;
    d_fe_family = LAGRANGE;
    d_fe_order = INVALID_ORDER;
    d_do_log = false;

    // Initialize function pointers to NULL.
    d_coordinate_mapping_fcns.resize(d_num_parts,NULL);
    d_coordinate_mapping_fcn_ctxs.resize(d_num_parts,NULL);
    d_PK1_stress_fcns.resize(d_num_parts,NULL);
    d_PK1_stress_fcn_systems.resize(d_num_parts);
    d_PK1_stress_fcn_ctxs.resize(d_num_parts,NULL);

    // Determine whether we should use first-order or second-order shape
    // functions for each part of the structure.
    bool mesh_has_first_order_elems  = false;
    bool mesh_has_second_order_elems = false;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const MeshBase& mesh = *meshes[part];
        MeshBase::const_element_iterator       el_it  = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for ( ; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            mesh_has_first_order_elems  = mesh_has_first_order_elems  || elem->default_order() == FIRST ;
            mesh_has_second_order_elems = mesh_has_second_order_elems || elem->default_order() == SECOND;
        }
    }
    mesh_has_first_order_elems  = SAMRAI_MPI::maxReduction(mesh_has_first_order_elems );
    mesh_has_second_order_elems = SAMRAI_MPI::maxReduction(mesh_has_second_order_elems);
    if (( mesh_has_first_order_elems &&  mesh_has_second_order_elems) ||
        (!mesh_has_first_order_elems && !mesh_has_second_order_elems))
    {
        TBOX_ERROR(d_object_name << "::SimplifiedIBFEMethod():\n"
                   << "  all parts of FE mesh must contain only FIRST order elements or only SECOND order elements" << std::endl);
    }
    if (mesh_has_first_order_elems)
    {
        d_fe_order = FIRST;
    }
    if (mesh_has_second_order_elems)
    {
        d_fe_order = SECOND;
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Report configuration.
    pout << "\n";
    pout << d_object_name << ": using " << Utility::enum_to_string<Order>(d_fe_order) << " order " << Utility::enum_to_string<FEFamily>(d_fe_family) << " finite elements.\n";
    pout << "\n";

    // Create the FE data managers that manage mappings between the FE mesh
    // parts and the Cartesian grid.
    d_meshes = meshes;
    d_equation_systems.resize(d_num_parts, NULL);
    d_fe_data_managers.resize(d_num_parts, NULL);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Create FE data managers.
        std::ostringstream manager_stream;
        manager_stream << "SimplifiedIBFEMethod FEDataManager::" << part;
        const std::string& manager_name = manager_stream.str();
        d_fe_data_managers[part] = FEDataManager::getManager(manager_name, "IB_4", "IB_4", d_use_consistent_mass_matrix);
        d_ghosts = IntVector<NDIM>::max(d_ghosts,d_fe_data_managers[part]->getGhostCellWidth());

        // Create FE equation systems object.
        d_equation_systems[part] = new EquationSystems(*d_meshes[part]);
        EquationSystems* equation_systems = d_equation_systems[part];
        d_fe_data_managers[part]->setEquationSystems(equation_systems, max_level_number-1);

        // Create FE systems and corresponding variables.
        d_fe_data_managers[part]->COORDINATES_SYSTEM_NAME = COORDS_SYSTEM_NAME;
        System& X_system = equation_systems->add_system<System>(COORDS_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "X_" << d;
            X_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        System& X_mapping_system = equation_systems->add_system<System>(COORD_MAPPING_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "dX_" << d;
            X_mapping_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        System& U_system = equation_systems->add_system<System>(VELOCITY_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "U_" << d;
            U_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }
    }

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time     = std::numeric_limits<double>::quiet_NaN();
    d_half_time    = std::numeric_limits<double>::quiet_NaN();

    d_is_initialized = false;
    return;
}// commonConstructor

void
SimplifiedIBFEMethod::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
    if (!is_from_restart)
    {
        if (db->isBool("use_consistent_mass_matrix")) d_use_consistent_mass_matrix = db->getBool("use_consistent_mass_matrix");
        if (db->isInteger("min_ghost_cell_width"))
        {
            d_ghosts = db->getInteger("min_ghost_cell_width");
        }
        else if (db->isDouble("min_ghost_cell_width"))
        {
            d_ghosts = static_cast<int>(std::ceil(db->getDouble("min_ghost_cell_width")));
        }
    }
    if      (db->keyExists("do_log"        )) d_do_log = db->getBool("do_log"        );
    else if (db->keyExists("enable_logging")) d_do_log = db->getBool("enable_logging");
    return;
}// getFromInput

void
SimplifiedIBFEMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("SIMPLIFIED_IBFE_METHOD_VERSION");
    if (ver != SIMPLIFIED_IBFE_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    db->getIntegerArray("d_ghosts", d_ghosts, NDIM);
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    d_fe_family = Utility::string_to_enum<FEFamily>(db->getString("d_fe_family"));
    d_fe_order = Utility::string_to_enum<Order>(db->getString("d_fe_order"));
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
