// Filename: CirculationModel.C
// Created on 20 Aug 2007 by Boyce Griffith

#include "CirculationModel.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <PatchLevel.h>
#include <tbox/SAMRAI_MPI.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CirculationModel::CirculationModel(
    const string& /*object_name*/,
    Pointer<Database> /*input_db*/)
{
    d_P_mean = 0.0;
    d_radius = 16.0*0.0416667;
    return;
}// CirculationModel

CirculationModel::~CirculationModel()
{
    return;
}// ~CirculationModel

void
CirculationModel::advanceTimeDependentData(
    const double data_time,
    const double dt,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int /*U_idx*/,
    const int P_idx,
    const int wgt_cc_idx,
    const int /*wgt_sc_idx*/)
{
    // Compute the mean pressure in the center of the domain.
    d_P_mean = 0.0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM,double> >      P_data = patch->getPatchData(     P_idx);
            Pointer<CellData<NDIM,double> > wgt_cc_data = patch->getPatchData(wgt_cc_idx);
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const double* const x_lower = pgeom->getXLower();
            const double* const dx = pgeom->getDx();
            double X[NDIM];
            for (Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const Index<NDIM>& i = b();
                double r_sq = 0.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d]*(double(i(d)-patch_box.lower(d))+0.5);
                    r_sq += pow(X[d] - 0.0, 2.0);
                }
                const double r = sqrt(r_sq);
                if (r <= d_radius) d_P_mean += (*P_data)(i)*(*wgt_cc_data)(i)/(M_PI*d_radius*d_radius);
            }
        }
    }
    SAMRAI_MPI::sumReduction(&d_P_mean,1);
    // R Q + L dQ/dt = P_rsvr - P.
    double R_src = 79.993432*1.5;
    double L_src = 1.0e0;
    double P_rsvr = 700.0*1333.2239; // 1000.0*1333.2239; //(100.0+900.0*0.5*(-cos(2.0*M_PI*data_time/0.01)+1.0))*1333.2239;
    d_Q = (2.0*P_rsvr*dt-2.0*d_P_mean*dt-R_src*d_Q*dt+2.0*L_src*d_Q)/(R_src*dt+2.0*L_src);
    return;
}// advanceTimeDependentData

bool
CirculationModel::isTimeDependent() const
{
    return true;
}// isTimeDependent

void
CirculationModel::setDataOnPatch(
    int data_idx,
    Pointer<Variable<NDIM> > /*var*/,
    Pointer<Patch<NDIM> > patch,
    double /*data_time*/,
    bool /*initial_time*/,
    Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    Pointer<CellData<NDIM,double> > Q_data = patch->getPatchData(data_idx);
    Q_data->fillAll(0.0);
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const Box<NDIM>& patch_box = patch->getBox();
    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();
    double X[NDIM];
    for (Box<NDIM>::Iterator b(patch_box); b; b++)
    {
        const Index<NDIM>& i = b();
        double r_sq = 0.0;
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = x_lower[d] + dx[d]*(double(i(d)-patch_box.lower(d))+0.5);
            r_sq += pow(X[d] - 0.0, 2.0);
        }
        const double r = sqrt(r_sq);
        if (r <= d_radius) (*Q_data)(i) = d_Q/(M_PI*d_radius*d_radius);
    }
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
