// Filename: WENOMOLOpsAndSoln.C
// Last modified: <01.Apr.2008 17:11:36 griffith@box221.cims.nyu.edu>
// Created on 04 Jan 2008 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "WENOMOLOps.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/PatchMathOps.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <FaceData.h>
#include <Index.h>
#include <MethodOfLinesIntegrator.h>

// C++ STDLIB INCLUDES
#include <algorithm>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define WENO_CONVECTIVE_FLUXES_F77 F77_FUNC_(weno_convective_fluxes2d, WENO_CONVECTIVE_FLUXES2D)
#endif

extern "C"
{
    void
    WENO_CONVECTIVE_FLUXES_F77(
        double* fface0, double* fface1, const int& fface_gcw,
        double* ffwrd0, double* ffwrd1, const int& ffwrd_gcw,
        double* frevr0, double* frevr1, const int& frevr_gcw,
        double* fplus0, double* fplus1, const int& fplus_gcw,
        double* fminus0, double* fminus1, const int& fminus_gcw,
        double* F, const int& F_gcw,
        const double* Q, const int& Q_gcw,
        const double* U, const int& U_gcw,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1);
}

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int GHOSTS = 4;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

WENOMOLOps::WENOMOLOps(
    const string& object_name,
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
    tbox::Pointer<tbox::Database> input_db)
    : d_object_name(object_name),
      d_visit_writer(NULL),
      d_grid_geom(grid_geom),
      d_X(NDIM),
      d_init_type("GAUSSIAN"),
      d_gaussian_kappa(0.01),
      d_disk_r(0.15),
      d_zalesak_r(0.15),
      d_zalesak_slot_w(0.025),
      d_zalesak_slot_l(0.1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!grid_geom.isNull());
#endif

    // Default initial values.
    const double* const XUpper = d_grid_geom->getXUpper();
    const double* const XLower = d_grid_geom->getXLower();
    for (int d = 0; d < NDIM; ++d)
    {
        d_X[d] = XLower[d] + 0.5*(XUpper[d] - XLower[d]);
    }

    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Initialize variables.
    d_Q_var = new pdat::CellVariable<NDIM,double>("Q",1);
    d_U_var = new pdat::CellVariable<NDIM,double>("U",NDIM);
    d_F_var = new pdat::CellVariable<NDIM,double>("F",1);
    d_flux_var = new pdat::FaceVariable<NDIM,double>("flux",1);
    return;
}// WENOMOLOps

WENOMOLOps::~WENOMOLOps()
{
    // intentionally blank
    return;
}// ~WENOMOLOps

void
WENOMOLOps::registerVisItDataWriter(
    tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_writer)
{
    d_visit_writer = visit_writer;
    return;
}// registerVisItDataWriter

void
WENOMOLOps::registerModelVariables(
    algs::MethodOfLinesIntegrator<NDIM>* integrator)
{
    // Register variables with the integrator.
    integrator->registerVariable(d_Q_var,GHOSTS,
                                 algs::MethodOfLinesIntegrator<NDIM>::SOLN,
                                 d_grid_geom,
                                 "CONSERVATIVE_COARSEN",
                                 "LINEAR_REFINE");

    integrator->registerVariable(d_U_var,GHOSTS,
                                 algs::MethodOfLinesIntegrator<NDIM>::SOLN,
                                 d_grid_geom,
                                 "CONSERVATIVE_COARSEN",
                                 "LINEAR_REFINE");

    integrator->registerVariable(d_F_var,GHOSTS,
                                 algs::MethodOfLinesIntegrator<NDIM>::RHS,
                                 d_grid_geom,
                                 "NO_COARSEN",
                                 "NO_REFINE");

    integrator->registerVariable(d_flux_var,1,
                                 algs::MethodOfLinesIntegrator<NDIM>::RHS,
                                 d_grid_geom,
                                 "NO_COARSEN",
                                 "NO_REFINE");

    // Register variables with the VisIt data writer.
    if (!(d_visit_writer.isNull()))
    {
        hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
        int Q_idx = var_db->mapVariableAndContextToIndex(d_Q_var, getInteriorContext());
        d_visit_writer->registerPlotQuantity(d_Q_var->getName(), "SCALAR", Q_idx);
    }
    return;
}// registerModelVariables

void
WENOMOLOps::initializeDataOnPatch(
        hier::Patch<NDIM>& patch,
        const double time,
        const bool initial_time) const

{
    hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
    int Q_idx = var_db->mapVariableAndContextToIndex(d_Q_var, getInteriorContext());
    tbox::Pointer< pdat::CellData<NDIM,double> > Q_data = patch.getPatchData(Q_idx);
    int U_idx = var_db->mapVariableAndContextToIndex(d_U_var, getInteriorContext());
    tbox::Pointer< pdat::CellData<NDIM,double> > U_data = patch.getPatchData(U_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_data.isNull());
    TBOX_ASSERT(!U_data.isNull());
#endif
    U_data->fill(-2.0,0);
    U_data->fill(+1.0,1);

    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    double r_squared;
    double X[NDIM];

    Q_data->fillAll(0.0);

    if (d_init_type == "GAUSSIAN")
    {
        for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            // NOTE: This assumes the lattice of Gaussians is being advected in
            // the unit square.
            int offset[NDIM];
            for (offset[0] = -2; offset[0] <= 2; ++(offset[0]))
            {
#if (NDIM>1)
                for (offset[1] = -2; offset[1] <= 2; ++(offset[1]))
                {
#endif
#if (NDIM>2)
                    for (offset[2] = -2; offset[2] <= 2; ++(offset[2]))
                    {
#endif
                        r_squared = 0.0;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X[d] = XLower[d] +
                                dx[d]*(double(i(d)-patch_lower(d))+0.5);
                            r_squared += pow(
                                X[d]-(d_X[d]+double(offset[d])),2.0);
                        }

                        (*Q_data)(i) +=
                            exp(-r_squared/(4.0*d_gaussian_kappa))/
                            pow(4.0*M_PI*d_gaussian_kappa,
                                0.5*double(NDIM));
#if (NDIM>2)
                    }
#endif
#if (NDIM>1)
                }
#endif
            }
        }
    }
    else if (d_init_type == "DISK")
    {
        for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            r_squared = 0.0;
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = XLower[d] +
                    dx[d]*(double(i(d)-patch_lower(d))+0.5);
                r_squared += pow((X[d]-d_X[d]),2.0);
            }
            if (sqrt(r_squared) > d_disk_r)
            {
                (*Q_data)(i) = 0.0;
            }
            else
            {
                (*Q_data)(i) = 1.0;
            }
        }
    }
    else if (d_init_type == "ZALESAK")
    {
        for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            r_squared = 0.0;
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = XLower[d] +
                    dx[d]*(double(i(d)-patch_lower(d))+0.5);
                r_squared += pow((X[d]-d_X[d]),2.0);
            }
            if ((sqrt(r_squared) > d_zalesak_r) ||
                ((abs(X[0] - d_X[0]) < d_zalesak_slot_w) &&
                 (X[1] - d_X[1]) < d_zalesak_slot_l))
            {
                (*Q_data)(i) = 0.0;
            }
            else
            {
                (*Q_data)(i) = 1.0;
            }
        }
    }
    else if (d_init_type == "SINUSOIDAL")
    {
        for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = XLower[d] +
                    dx[d]*(double(i(d)-patch_lower(d))+0.5);
            }
            (*Q_data)(i) = sin(X[0]);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeDataOnPatch()\n"
                   << "  invalid initialization type " << d_init_type << "\n");
    }
    return;
}// initializeDataOnPatch

double
WENOMOLOps::computeStableDtOnPatch(
    hier::Patch<NDIM>& patch,
    const double time) const
{
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();

    static const double U_max = 2.0;
    double dt = numeric_limits<double>::max();
    for (int d = 0; d < NDIM; ++d)
    {
        dt = min(dt, dx[d]/U_max);
    }

    return 0.475*dt;
}// computeStableDtOnPatch

void
WENOMOLOps::singleStep(
    hier::Patch<NDIM>& patch,
    const double dt,
    const double alpha_1,
    const double alpha_2,
    const double beta) const
{
    tbox::Pointer<pdat::CellData<NDIM,double> > Q_updated =
        patch.getPatchData(d_Q_var, getInteriorWithGhostsContext());
    const int Q_gcw = Q_updated->getGhostCellWidth().max();

    tbox::Pointer<pdat::CellData<NDIM,double> > Q_fixed =
        patch.getPatchData(d_Q_var, getInteriorContext());

    tbox::Pointer<pdat::CellData<NDIM,double> > U_updated =
        patch.getPatchData(d_U_var, getInteriorWithGhostsContext());
    const int U_gcw = U_updated->getGhostCellWidth().max();

    tbox::Pointer<pdat::CellData<NDIM,double> > U_fixed =
        patch.getPatchData(d_U_var, getInteriorContext());

    tbox::Pointer<pdat::CellData<NDIM,double> > F =
        patch.getPatchData(d_F_var, getInteriorContext());
    const int F_gcw = F->getGhostCellWidth().max();

    tbox::Pointer<pdat::FaceData<NDIM,double> > flux =
        patch.getPatchData(d_flux_var, getInteriorContext());
    const int flux_gcw = flux->getGhostCellWidth().max();

    tbox::Pointer<pdat::FaceData<NDIM,double> > flux_fwrd =
        new pdat::FaceData<NDIM,double>(
            flux->getBox(), flux->getDepth(), flux->getGhostCellWidth());
    const int flux_fwrd_gcw = flux_fwrd->getGhostCellWidth().max();

    tbox::Pointer<pdat::FaceData<NDIM,double> > flux_revr =
        new pdat::FaceData<NDIM,double>(
            flux->getBox(), flux->getDepth(), flux->getGhostCellWidth());
    const int flux_revr_gcw = flux_revr->getGhostCellWidth().max();

    tbox::Pointer<pdat::FaceData<NDIM,double> > flux_plus =
        new pdat::FaceData<NDIM,double>(
            flux->getBox(), flux->getDepth(), flux->getGhostCellWidth());
    const int flux_plus_gcw = flux_plus->getGhostCellWidth().max();

    tbox::Pointer<pdat::FaceData<NDIM,double> > flux_minus =
        new pdat::FaceData<NDIM,double>(
            flux->getBox(), flux->getDepth(), flux->getGhostCellWidth());
    const int flux_minus_gcw = flux_minus->getGhostCellWidth().max();

    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    const hier::Index<NDIM>& patch_upper = patch_box.upper();

    WENO_CONVECTIVE_FLUXES_F77(
        flux->getPointer(0), flux->getPointer(1), flux_gcw,
        flux_fwrd->getPointer(0), flux_fwrd->getPointer(1), flux_fwrd_gcw,
        flux_revr->getPointer(0), flux_revr->getPointer(1), flux_revr_gcw,
        flux_plus->getPointer(0), flux_plus->getPointer(1), flux_plus_gcw,
        flux_minus->getPointer(0), flux_minus->getPointer(1), flux_minus_gcw,
        F->getPointer(), F_gcw,
        Q_updated->getPointer(), Q_gcw,
        U_updated->getPointer(), U_gcw,
        patch_lower(0), patch_upper(0),
        patch_lower(1), patch_upper(1));

    IBTK::PatchMathOps patch_math_ops;
    patch_math_ops.div(F,
                       -1.0, flux,
                       0.0, tbox::Pointer<pdat::CellData<NDIM,double> >(NULL),
                       tbox::Pointer<SAMRAI::hier::Patch<NDIM> >(&patch,false));

    for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
    {
        const hier::Index<NDIM>& i = ic();
        (*Q_updated)(i) = alpha_1*(*Q_fixed)(i) + alpha_2*(*Q_updated)(i) + beta*dt*(*F)(i);
    }
    return;
}// singleStep

void
WENOMOLOps::setPhysicalBoundaryConditions(
    hier::Patch<NDIM>& patch,
    const double fill_time,
    const hier::IntVector<NDIM>& ghost_width_to_fill)
{
    TBOX_ERROR("physical boundary conditions are not implemented\n");
    return;
}// setPhysicalBoundaryConditions

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
WENOMOLOps::getFromInput(
    tbox::Pointer<tbox::Database> db)
{
    if (!db.isNull())
    {
        if (db->keyExists("X"))
        {
            d_X = db->getDoubleArray("X");
        }

        d_init_type = db->getStringWithDefault("init_type",d_init_type);

        if (d_init_type == "GAUSSIAN")
        {
            d_gaussian_kappa = db->
                getDoubleWithDefault("kappa",d_gaussian_kappa);
        }
        else if (d_init_type == "DISK")
        {
            d_disk_r = db->
                getDoubleWithDefault("disk_r",d_disk_r);
        }
        else if (d_init_type == "ZALESAK")
        {
            d_zalesak_r = db->
                getDoubleWithDefault("zalesak_r",d_zalesak_r);
            d_zalesak_slot_w = db->
                getDoubleWithDefault("zalesak_slot_w",d_zalesak_slot_w);
            d_zalesak_slot_l = db->
                getDoubleWithDefault("zalesak_slot_l",d_zalesak_slot_l);
        }
        else if (d_init_type == "SINUSOIDAL")
        {
            // intentionally blank
        }
        else
        {
            TBOX_ERROR(d_object_name << "::getFromInput()\n"
                       << "  invalid initialization type " << d_init_type << "\n");
        }
    }
    return;
}// getFromInput

//////////////////////////////////////////////////////////////////////////////
