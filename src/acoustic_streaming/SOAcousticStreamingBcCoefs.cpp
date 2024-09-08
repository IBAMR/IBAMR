// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/SOAcousticStreamingBcCoefs.h"

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "SideData.h"
#include "SideIndex.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{

inline hier::Index<NDIM>
get_shift(int dir, int shift)
{
    SAMRAI::hier::Index<NDIM> iv(0);
    iv(dir) = shift;
    return iv;
} // get_shift

static const int EXTENSIONS_FILLABLE = 128;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
SOAcousticStreamingBcCoefs::setSOVelocityComponent(int U2_comp)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(U2_comp >= 0 && U2_comp < NDIM);
#endif

    d_U2_comp = U2_comp;
    return;
} // setSOVelocityComponent

void
SOAcousticStreamingBcCoefs::setFOVelocityPressurePatchDataIndices(int U1_real_idx,
                                                                  int U1_imag_idx,
                                                                  int p1_real_idx,
                                                                  int p1_imag_idx)
{
    d_U1_real_idx = U1_real_idx;
    d_U1_imag_idx = U1_imag_idx;
    d_p1_real_idx = p1_real_idx;
    d_p1_imag_idx = p1_imag_idx;
    return;
} // setFOVelocityPressurePatchDataIndices

void
SOAcousticStreamingBcCoefs::setDensityPatchDataIndex(int rho_idx)
{
    d_rho_idx = rho_idx;
    return;
} // setDensityPatchDataIndex

void
SOAcousticStreamingBcCoefs::setSoundSpeed(double sound_speed)
{
    d_sound_speed = sound_speed;
    return;
} // setSoundSpeed

void
SOAcousticStreamingBcCoefs::setAcousticAngularFrequency(double omega)
{
    d_acoustic_freq = omega;
    return;
} // setAcousticAngularFrequency

void
SOAcousticStreamingBcCoefs::useStokesDriftVelocityForm(bool use_stokes_drift_vel)
{
    d_use_stokes_drift_velocity = use_stokes_drift_vel;
    return;
} // useStokesDriftVelocityForm

void
SOAcousticStreamingBcCoefs::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                       Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                       Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                       const Pointer<Variable<NDIM> >& /*variable*/,
                                       const Patch<NDIM>& patch,
                                       const BoundaryBox<NDIM>& bdry_box,
                                       double fill_time) const
{
    if (d_use_stokes_drift_velocity)
    {
        setBcCoefsFromStokesDrift(acoef_data, bcoef_data, gcoef_data, patch, bdry_box, fill_time);
    }
    else
    {
        setBcCoefsFromFOMassFlux(acoef_data, bcoef_data, gcoef_data, patch, bdry_box, fill_time);
    }

    return;

} // setBcCoefs

IntVector<NDIM>
SOAcousticStreamingBcCoefs::numberOfExtensionsFillable() const
{
    return EXTENSIONS_FILLABLE;
} // numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SOAcousticStreamingBcCoefs::setBcCoefsFromFOMassFlux(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                                     Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                                     Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                                     const Patch<NDIM>& patch,
                                                     const BoundaryBox<NDIM>& bdry_box,
                                                     double fill_time) const
{
    // Patch data for first-order velocity and density
    Pointer<SideData<NDIM, double> > U1_real_data = patch.getPatchData(d_U1_real_idx);
    Pointer<SideData<NDIM, double> > U1_imag_data = patch.getPatchData(d_U1_imag_idx);
    Pointer<CellData<NDIM, double> > p1_real_data = patch.getPatchData(d_p1_real_idx);
    Pointer<CellData<NDIM, double> > p1_imag_data = patch.getPatchData(d_p1_imag_idx);
    Pointer<SideData<NDIM, double> > rho_data = patch.getPatchData(d_rho_idx);

    // Loop over the boundary box and set the coefficients.
    const unsigned int location_index = bdry_box.getLocationIndex();
    const int bdry_normal_axis = location_index / 2;
    const bool normal_comp = (d_U2_comp == bdry_normal_axis);

    const Box<NDIM>& bc_coef_box = (acoef_data ? acoef_data->getBox() :
                                    bcoef_data ? bcoef_data->getBox() :
                                    gcoef_data ? gcoef_data->getBox() :
                                                 Box<NDIM>());
#if !defined(NDEBUG)
    TBOX_ASSERT(!acoef_data || bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(!bcoef_data || bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(!gcoef_data || bc_coef_box == gcoef_data->getBox());
#endif

    for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
    {
        const hier::Index<NDIM>& i = b();

        if (acoef_data) (*acoef_data)(i, 0) = 1.0;
        if (bcoef_data) (*bcoef_data)(i, 0) = 0.0;

        //     u2 = -1/rho0 <rho1 u1>
        // ==> u2 = -1/rho0 * 1/2 *(rho1,real u1,real + rho1,imag u1,imag)
        if (gcoef_data)
        {
            if (normal_comp)
            {
                const SideIndex<NDIM> is(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                const hier::Index<NDIM> shift_w = get_shift(bdry_normal_axis, -1);

                const double& rho0 = (*rho_data)(is);
                const double& U1_real = (*U1_real_data)(is);
                const double& U1_imag = (*U1_imag_data)(is);
                const double rho1_real =
                    0.5 * ((*p1_real_data)(i) + (*p1_real_data)(i + shift_w)) / (d_sound_speed * d_sound_speed);
                const double rho1_imag =
                    0.5 * ((*p1_imag_data)(i) + (*p1_imag_data)(i + shift_w)) / (d_sound_speed * d_sound_speed);

                (*gcoef_data)(i, 0) = (-1.0 / rho0) * 0.5 * (rho1_real * U1_real + rho1_imag * U1_imag);
            }
            else
            {
                const hier::Index<NDIM> shift_w = get_shift(bdry_normal_axis, -1);
                const hier::Index<NDIM> shift_s = get_shift(d_U2_comp, -1);
                const hier::Index<NDIM> shift_sw = shift_w + shift_s;

                const SideIndex<NDIM> is(i, d_U2_comp, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> is_w(i + shift_w, d_U2_comp, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> is_n(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> is_s(i + shift_s, bdry_normal_axis, SideIndex<NDIM>::Lower);

                const double U1_real = 0.5 * ((*U1_real_data)(is) + (*U1_real_data)(is_w));
                const double U1_imag = 0.5 * ((*U1_imag_data)(is) + (*U1_imag_data)(is_w));
                const double rho0 = 0.5 * ((*rho_data)(is) + (*rho_data)(is_w));
                const double rho1_real = 0.25 *
                                         ((*p1_real_data)(i) + (*p1_real_data)(i + shift_w) +
                                          (*p1_real_data)(i + shift_s) + (*p1_real_data)(i + shift_sw)) /
                                         (d_sound_speed * d_sound_speed);
                const double rho1_imag = 0.25 *
                                         ((*p1_imag_data)(i) + (*p1_imag_data)(i + shift_w) +
                                          (*p1_imag_data)(i + shift_s) + (*p1_imag_data)(i + shift_sw)) /
                                         (d_sound_speed * d_sound_speed);

                (*gcoef_data)(i, 0) = (-1.0 / rho0) * 0.5 * (rho1_real * U1_real + rho1_imag * U1_imag);
            }
        }
    }
    return;
} // setBcCoefsFromFOMassFlux

void
SOAcousticStreamingBcCoefs::setBcCoefsFromStokesDrift(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                                      Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                                      Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                                      const Patch<NDIM>& patch,
                                                      const BoundaryBox<NDIM>& bdry_box,
                                                      double fill_time) const
{
    // Patch data for first-order velocity and density
    Pointer<SideData<NDIM, double> > Ur = patch.getPatchData(d_U1_real_idx);
    Pointer<SideData<NDIM, double> > Ui = patch.getPatchData(d_U1_imag_idx);
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Loop over the boundary box and set the coefficients.
    const unsigned int location_index = bdry_box.getLocationIndex();
    const int bdry_normal_axis = location_index / 2;
    const bool normal_comp = (d_U2_comp == bdry_normal_axis);

    const Box<NDIM>& bc_coef_box = (acoef_data ? acoef_data->getBox() :
                                    bcoef_data ? bcoef_data->getBox() :
                                    gcoef_data ? gcoef_data->getBox() :
                                                 Box<NDIM>());
#if !defined(NDEBUG)
    TBOX_ASSERT(!acoef_data || bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(!bcoef_data || bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(!gcoef_data || bc_coef_box == gcoef_data->getBox());
#endif

    for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
    {
        const hier::Index<NDIM>& i = b();

        if (acoef_data) (*acoef_data)(i, 0) = 1.0;
        if (bcoef_data) (*bcoef_data)(i, 0) = 0.0;

        //     u2 = -u_sd
        if (gcoef_data)
        {
            if (normal_comp)
            {
                const SideIndex<NDIM> is(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> is_e(i, bdry_normal_axis, SideIndex<NDIM>::Upper);
                const SideIndex<NDIM> is_w(
                    i + get_shift(bdry_normal_axis, -1), bdry_normal_axis, SideIndex<NDIM>::Lower);
                const double g_normal =
                    ((*Ui)(is) * ((*Ur)(is_e) - (*Ur)(is_w)) - (*Ur)(is) * ((*Ui)(is_e) - (*Ui)(is_w))) /
                    (2.0 * dx[bdry_normal_axis]);

                double g_tangential = 0.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    if (d == bdry_normal_axis) continue;

                    const SideIndex<NDIM> is_n(i + get_shift(d, 1), bdry_normal_axis, SideIndex<NDIM>::Lower);
                    const SideIndex<NDIM> is_s(i + get_shift(d, -1), bdry_normal_axis, SideIndex<NDIM>::Lower);

                    const SideIndex<NDIM> is_e(i, d, SideIndex<NDIM>::Lower);
                    const SideIndex<NDIM> is_ne(i, d, SideIndex<NDIM>::Upper);
                    const SideIndex<NDIM> is_w(i + get_shift(bdry_normal_axis, -1), d, SideIndex<NDIM>::Lower);
                    const SideIndex<NDIM> is_nw(i + get_shift(bdry_normal_axis, -1), d, SideIndex<NDIM>::Upper);

                    g_tangential += 0.25 * ((*Ui)(is_e) + (*Ui)(is_ne) + (*Ui)(is_w) + (*Ui)(is_nw)) *
                                    ((*Ur)(is_n) - (*Ur)(is_s)) / (2 * dx[d]);
                    g_tangential -= 0.25 * ((*Ur)(is_e) + (*Ur)(is_ne) + (*Ur)(is_w) + (*Ur)(is_nw)) *
                                    ((*Ui)(is_n) - (*Ui)(is_s)) / (2 * dx[d]);
                }

                (*gcoef_data)(i, 0) = -(g_normal + g_tangential) / (2.0 * d_acoustic_freq);
            }
            else
            {
                const SideIndex<NDIM> is_e(i, d_U2_comp, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> is_w(i + get_shift(bdry_normal_axis, -1), d_U2_comp, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> is_n(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> is_s(i + get_shift(d_U2_comp, -1), bdry_normal_axis, SideIndex<NDIM>::Lower);

                const double g_normal =
                    0.5 * ((*Ui)(is_n) + (*Ui)(is_s)) * (((*Ur)(is_e) - (*Ur)(is_w)) / dx[bdry_normal_axis]) -
                    0.5 * ((*Ur)(is_n) + (*Ur)(is_s)) * (((*Ui)(is_e) - (*Ui)(is_w)) / dx[bdry_normal_axis]);

                double g_tangential = 0.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    if (d == bdry_normal_axis) continue;

                    const SideIndex<NDIM> is_ne(i, d_U2_comp, SideIndex<NDIM>::Upper);
                    const SideIndex<NDIM> is_se(i + get_shift(d, -1), d_U2_comp, SideIndex<NDIM>::Lower);
                    const SideIndex<NDIM> is_nw(i + get_shift(bdry_normal_axis, -1), d_U2_comp, SideIndex<NDIM>::Upper);
                    const SideIndex<NDIM> is_sw(
                        i + get_shift(bdry_normal_axis, -1) + get_shift(d, -1), d_U2_comp, SideIndex<NDIM>::Lower);

                    g_tangential += 0.5 * ((*Ui)(is_e) + (*Ui)(is_w)) * 0.5 *
                                    ((*Ur)(is_ne) - (*Ur)(is_se) + (*Ur)(is_nw) - (*Ur)(is_sw)) / (2 * dx[d]);
                    g_tangential -= 0.5 * ((*Ur)(is_e) + (*Ur)(is_w)) * 0.5 *
                                    ((*Ui)(is_ne) - (*Ui)(is_se) + (*Ui)(is_nw) - (*Ui)(is_sw)) / (2 * dx[d]);
                }

                (*gcoef_data)(i, 0) = -(g_normal + g_tangential) / (2.0 * d_acoustic_freq);
            }
        }
    }
    return;
} // setBcCoefsFromStokesDrift

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
