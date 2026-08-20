#include <ibamr/RMTMethod.h>

#include <ibtk/IBTK_MPI.h>

#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <Index.h>
#include <Patch.h>
#include <PatchLevel.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <ibamr/RMTMethod.h>
#include <ibamr/RMTForceFunction.h>

RMTMethod::RMTMethod(const std::string& object_name,
                     Pointer<PatchHierarchy<NDIM> > hierarchy,
                     Pointer<CellVariable<NDIM, double> > phi_var,
                     Pointer<CellVariable<NDIM, double> > xi0_var,
                     Pointer<CellVariable<NDIM, double> > xi1_var,
                     int phi_idx,
                     int xi0_idx,
                     int xi1_idx,
                     double band_width,
                     int num_passes,
                     int reinit_num_iter,
                     double reinit_dtau,
                     double GGS,
                     double nu_s,
                     double rho_s,
                     double g0y,
                     double Tramp)
    : d_object_name(object_name),
      d_hierarchy(hierarchy),
      d_phi_idx(phi_idx),
      d_xi0_idx(xi0_idx),
      d_xi1_idx(xi1_idx),
      d_band_width(band_width),
      d_num_passes(num_passes),
      d_reinit_num_iter(reinit_num_iter),
      d_reinit_dtau(reinit_dtau)
{

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

d_phi0_var = new CellVariable<NDIM, double>("PHI0");
d_w_var = new CellVariable<NDIM, double>("W_FROM_XI");
d_phiw_err_var = new CellVariable<NDIM, double>("PHI_MINUS_W");

d_phi0_idx = var_db->registerVariableAndContext(
    d_phi0_var,
    var_db->getContext("PHI0_CONTEXT"),
    IntVector<NDIM>(1));

d_w_idx = var_db->registerVariableAndContext(
    d_w_var,
    var_db->getContext("W_FROM_XI_CONTEXT"),
    IntVector<NDIM>(1));

d_phiw_err_idx = var_db->registerVariableAndContext(
    d_phiw_err_var,
    var_db->getContext("PHI_MINUS_W_CTX"),
    IntVector<NDIM>(0));

    d_force_fcn =
        new RMTForceFunction("rmt_force",
                             phi_var,
                             xi0_var,
                             xi1_var,
                             GGS,
                             nu_s,
                             rho_s,
                             g0y,
                             Tramp);

    d_force_fcn->setPhi0AndWIndices(d_phi0_idx, d_w_idx);
}


Pointer<IBTK::CartGridFunction>
RMTMethod::getForceFunction() const
{
    return d_force_fcn;
}



void
RMTMethod::initializeLevelData(double init_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            if (!patch->checkAllocated(d_phi0_idx))
                patch->allocatePatchData(d_phi0_idx, init_time);

            if (!patch->checkAllocated(d_w_idx))
                patch->allocatePatchData(d_w_idx, init_time);

            if (!patch->checkAllocated(d_phiw_err_idx))
                patch->allocatePatchData(d_phiw_err_idx, init_time);

            auto phi0_data = Pointer<CellData<NDIM, double> >(patch->getPatchData(d_phi0_idx));
            auto phi_data  = Pointer<CellData<NDIM, double> >(patch->getPatchData(d_phi_idx));
            auto w_data    = Pointer<CellData<NDIM, double> >(patch->getPatchData(d_w_idx));

            TBOX_ASSERT(!phi0_data.isNull());
            TBOX_ASSERT(!phi_data.isNull());
            TBOX_ASSERT(!w_data.isNull());

            phi0_data->copy(*phi_data);
            w_data->fillAll(0.0);
        }
    }
}


void
RMTMethod::registerVisItDataWriter(Pointer<VisItDataWriter<NDIM> > visit_data_writer)
{
    if (!visit_data_writer.isNull())
    {
        visit_data_writer->registerPlotQuantity("W_FROM_XI", "SCALAR", d_w_idx, 0, 1, "CELL");
        visit_data_writer->registerPlotQuantity("PHI_MINUS_W", "SCALAR", d_phiw_err_idx, 0, 1, "CELL");
    }
}

void
RMTMethod::initializeForceFunction(int phi_cur,
                                   int phi_new,
                                   int phi_scr,
                                   int xi0_cur,
                                   int xi0_new,
                                   int xi0_scr,
                                   int xi1_cur,
                                   int xi1_new,
                                   int xi1_scr)
{
    d_force_fcn->setIndices(phi_cur,
                            phi_new,
                            phi_scr,
                            xi0_cur,
                            xi0_new,
                            xi0_scr,
                            xi1_cur,
                            xi1_new,
                            xi1_scr);

    d_force_fcn->initializeGhostFill(d_hierarchy);
}



static inline double
clampDouble(const double x, const double a, const double b)
{
    return std::max(a, std::min(b, x));
}

static void
extrapolateXiNarrowBand(Pointer<PatchHierarchy<NDIM> > hierarchy,
                        int phi_idx,
                        int xi0_idx,
                        int xi1_idx,
                        double band_width,
                        int num_passes = -1)
{
    TBOX_ASSERT(NDIM == 2);

    const int finest_ln = hierarchy->getFinestLevelNumber();

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();

            auto phi = Pointer<CellData<NDIM, double> >(patch->getPatchData(phi_idx));
            auto xi0 = Pointer<CellData<NDIM, double> >(patch->getPatchData(xi0_idx));
            auto xi1 = Pointer<CellData<NDIM, double> >(patch->getPatchData(xi1_idx));

            TBOX_ASSERT(!phi.isNull());
            TBOX_ASSERT(!xi0.isNull());
            TBOX_ASSERT(!xi1.isNull());

            auto geom = Pointer<CartesianPatchGeometry<NDIM> >(patch->getPatchGeometry());
            TBOX_ASSERT(!geom.isNull());

            const double* dx = geom->getDx();
            const double* xlo = geom->getXLower();
            const hier::Index<NDIM>& ilo = box.lower();
            const hier::Index<NDIM>& ihi = box.upper();

            const int nx = ihi(0) - ilo(0) + 1;
            const int ny = ihi(1) - ilo(1) + 1;

            auto lid = [&](const CellIndex<NDIM>& ci) -> int { return (ci(0) - ilo(0)) + nx * (ci(1) - ilo(1)); };

            auto inside_box = [&](const CellIndex<NDIM>& ci) -> bool { return box.contains(ci); };

            auto cell_coord = [&](const CellIndex<NDIM>& ci, double& x, double& y)
            {
                x = xlo[0] + dx[0] * (double(ci(0) - ilo(0)) + 0.5);
                y = xlo[1] + dx[1] * (double(ci(1) - ilo(1)) + 0.5);
            };

            std::vector<int> known(nx * ny, 0);
            std::vector<int> target(nx * ny, 0);

            for (Box<NDIM>::Iterator it(box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double phiv = (*phi)(ci);

                if (phiv <= 0.0) known[lid(ci)] = 1;

                if (phiv <= band_width) target[lid(ci)] = 1;
            }

            const double r = std::sqrt(dx[0] * dx[0] + dx[1] * dx[1]);
            const double stencil_radius = 4.0 * r;
            const int si = std::max(2, int(std::ceil(stencil_radius / dx[0])));
            const int sj = std::max(2, int(std::ceil(stencil_radius / dx[1])));

            int max_passes = num_passes;
            if (max_passes <= 0) max_passes = int(std::ceil(band_width / std::min(dx[0], dx[1]))) + 2;

            auto solve3x3 = [](double A[3][3], double b[3], double c[3]) -> bool
            {
                double M[3][4];
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j) M[i][j] = A[i][j];
                    M[i][3] = b[i];
                }

                for (int k = 0; k < 3; ++k)
                {
                    int piv = k;
                    for (int i = k + 1; i < 3; ++i)
                        if (std::abs(M[i][k]) > std::abs(M[piv][k])) piv = i;

                    if (std::abs(M[piv][k]) < 1.0e-14) return false;

                    if (piv != k)
                    {
                        for (int j = k; j < 4; ++j) std::swap(M[k][j], M[piv][j]);
                    }

                    const double diag = M[k][k];
                    for (int j = k; j < 4; ++j) M[k][j] /= diag;

                    for (int i = 0; i < 3; ++i)
                    {
                        if (i == k) continue;
                        const double f = M[i][k];
                        for (int j = k; j < 4; ++j) M[i][j] -= f * M[k][j];
                    }
                }

                c[0] = M[0][3];
                c[1] = M[1][3];
                c[2] = M[2][3];
                return true;
            };

            for (int pass = 0; pass < max_passes; ++pass)
            {
                std::vector<int> new_known = known;
                std::vector<double> new_xi0(nx * ny, 0.0);
                std::vector<double> new_xi1(nx * ny, 0.0);
                int n_update = 0;

                for (Box<NDIM>::Iterator it(box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const int id = lid(ci);

                    if (!target[id] || known[id]) continue;

                    bool touches_known = false;
                    for (int dj = -1; dj <= 1; ++dj)
                    {
                        for (int di = -1; di <= 1; ++di)
                        {
                            if (di == 0 && dj == 0) continue;

                            hier::IntVector<NDIM> off(0);
                            off(0) = di;
                            off(1) = dj;
                            CellIndex<NDIM> nb = ci + off;

                            if (inside_box(nb) && known[lid(nb)]) touches_known = true;
                        }
                    }

                    if (!touches_known) continue;

                    double xc, yc;
                    cell_coord(ci, xc, yc);

                    double A[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
                    double b0[3] = { 0.0, 0.0, 0.0 };
                    double b1[3] = { 0.0, 0.0, 0.0 };

                    int npts = 0;

                    for (int jj = -sj; jj <= sj; ++jj)
                    {
                        for (int ii = -si; ii <= si; ++ii)
                        {
                            hier::IntVector<NDIM> off(0);
                            off(0) = ii;
                            off(1) = jj;
                            CellIndex<NDIM> nb = ci + off;

                            if (!inside_box(nb)) continue;
                            if (!known[lid(nb)]) continue;

                            double x, y;
                            cell_coord(nb, x, y);

                            const double dist = std::sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));

                            if (dist > stencil_radius) continue;

                            const double q[3] = { x, y, 1.0 };

                            for (int a = 0; a < 3; ++a)
                            {
                                for (int b = 0; b < 3; ++b) A[a][b] += q[a] * q[b];

                                b0[a] += q[a] * (*xi0)(nb);
                                b1[a] += q[a] * (*xi1)(nb);
                            }

                            ++npts;
                        }
                    }

                    if (npts < 3) continue;

                    double c0[3], c1[3];
                    double A0[3][3], A1[3][3];

                    for (int a = 0; a < 3; ++a)
                    {
                        for (int b = 0; b < 3; ++b)
                        {
                            A0[a][b] = A[a][b];
                            A1[a][b] = A[a][b];
                        }
                    }

                    if (!solve3x3(A0, b0, c0)) continue;
                    if (!solve3x3(A1, b1, c1)) continue;

                    new_xi0[id] = c0[0] * xc + c0[1] * yc + c0[2];
                    new_xi1[id] = c1[0] * xc + c1[1] * yc + c1[2];

                    new_known[id] = 1;
                    ++n_update;
                }

                for (Box<NDIM>::Iterator it(box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const int id = lid(ci);

                    if (!known[id] && new_known[id])
                    {
                        (*xi0)(ci) = new_xi0[id];
                        (*xi1)(ci) = new_xi1[id];
                    }
                }

                known.swap(new_known);

                if (n_update == 0) break;
            }

            //   delete phi_old;
            //     delete xi0_old;
            //     delete xi1_old;
        }
    }
}



static double
bilinearSampleCellDataPhysical(const SAMRAI::pdat::CellData<NDIM, double>& data,
                               const Box<NDIM>& box,
                               const double* x_lo,
                               const double* dx,
                               const hier::Index<NDIM>& ilo,
                               const double X,
                               const double Y)
{
    const double i_f = (X - x_lo[0]) / dx[0] + static_cast<double>(ilo(0)) - 0.5;
    const double j_f = (Y - x_lo[1]) / dx[1] + static_cast<double>(ilo(1)) - 0.5;

    const int i0 = static_cast<int>(std::floor(i_f));
    const int j0 = static_cast<int>(std::floor(j_f));
    const int i1 = i0 + 1;
    const int j1 = j0 + 1;

    CellIndex<NDIM> c00;
    c00(0) = i0;
    c00(1) = j0;
    CellIndex<NDIM> c10;
    c10(0) = i1;
    c10(1) = j0;
    CellIndex<NDIM> c01;
    c01(0) = i0;
    c01(1) = j1;
    CellIndex<NDIM> c11;
    c11(0) = i1;
    c11(1) = j1;

    CellIndex<NDIM> cc00;
    cc00(0) = std::max(box.lower(0), std::min(box.upper(0), c00(0)));
    cc00(1) = std::max(box.lower(1), std::min(box.upper(1), c00(1)));

    CellIndex<NDIM> cc10;
    cc10(0) = std::max(box.lower(0), std::min(box.upper(0), c10(0)));
    cc10(1) = std::max(box.lower(1), std::min(box.upper(1), c10(1)));

    CellIndex<NDIM> cc01;
    cc01(0) = std::max(box.lower(0), std::min(box.upper(0), c01(0)));
    cc01(1) = std::max(box.lower(1), std::min(box.upper(1), c01(1)));

    CellIndex<NDIM> cc11;
    cc11(0) = std::max(box.lower(0), std::min(box.upper(0), c11(0)));
    cc11(1) = std::max(box.lower(1), std::min(box.upper(1), c11(1)));

    const double tx = i_f - std::floor(i_f);
    const double ty = j_f - std::floor(j_f);

    const double v00 = data(cc00);
    const double v10 = data(cc10);
    const double v01 = data(cc01);
    const double v11 = data(cc11);

    const double v0 = (1.0 - tx) * v00 + tx * v10;
    const double v1 = (1.0 - tx) * v01 + tx * v11;

    return (1.0 - ty) * v0 + ty * v1;
}



static void
buildWFromXiUsingPhi0(Pointer<PatchHierarchy<NDIM> > hierarchy,
                      int w_idx,
                      int xi0_idx,
                      int xi1_idx,
                      int phi0_idx,
                      int phi_idx)
{
    TBOX_ASSERT(NDIM == 2);

    double max_err = 0.0;
    int imax = -1;
    int jmax = -1;
    double phi_at_max = 0.0;
    double w_at_max = 0.0;

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            auto w = Pointer<CellData<NDIM, double> >(patch->getPatchData(w_idx));
            auto xi0 = Pointer<CellData<NDIM, double> >(patch->getPatchData(xi0_idx));
            auto xi1 = Pointer<CellData<NDIM, double> >(patch->getPatchData(xi1_idx));
            auto phi0 = Pointer<CellData<NDIM, double> >(patch->getPatchData(phi0_idx));
            auto phi = Pointer<CellData<NDIM, double> >(patch->getPatchData(phi_idx));

            if (w.isNull() || xi0.isNull() || xi1.isNull() || phi0.isNull() || phi.isNull()) continue;

            const Box<NDIM>& box = patch->getBox();

            auto geom = Pointer<CartesianPatchGeometry<NDIM> >(patch->getPatchGeometry());
            const double* dx = geom->getDx();
            const double* x_lo = geom->getXLower();
            const hier::Index<NDIM>& ilo = box.lower();

            const double dx_max = std::max(dx[0], dx[1]);

            for (Box<NDIM>::Iterator it(box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                if ((*phi)(ci) > 2.0 * dx_max)
                {
                    (*w)(ci) = (*phi)(ci);
                    continue;
                }

                const double X = (*xi0)(ci);
                const double Y = (*xi1)(ci);

                const double wval = bilinearSampleCellDataPhysical(*phi0, box, x_lo, dx, ilo, X, Y);

                (*w)(ci) = wval;

                if (std::abs((*phi)(ci)) < 2.0 * dx_max)
                {
                    const double err = std::abs((*phi)(ci)-wval);

                    if (err > max_err)
                    {
                        max_err = err;
                        imax = ci(0);
                        jmax = ci(1);
                        phi_at_max = (*phi)(ci);
                        w_at_max = wval;
                    }
                }
            }
        }
    }

    pout << "[PHI-W MAX] "
         << " err=" << max_err << " i=" << imax << " j=" << jmax << " phi=" << phi_at_max << " w=" << w_at_max
         << std::endl;
}


static void
copyWToPhi(Pointer<PatchHierarchy<NDIM> > hierarchy, int phi_idx, int w_idx)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            auto phi = Pointer<CellData<NDIM, double> >(patch->getPatchData(phi_idx));

            auto w = Pointer<CellData<NDIM, double> >(patch->getPatchData(w_idx));

            if (phi.isNull() || w.isNull()) continue;

            phi->copy(*w);
        }
    }
}

static void
computePhiWError(Pointer<PatchHierarchy<NDIM> > hierarchy, int phi_idx, int w_idx)
{
    double err_max = 0.0;

    const int finest_ln = hierarchy->getFinestLevelNumber();

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            auto phi = Pointer<CellData<NDIM, double> >(patch->getPatchData(phi_idx));

            auto w = Pointer<CellData<NDIM, double> >(patch->getPatchData(w_idx));

            const Box<NDIM>& box = patch->getBox();

            for (Box<NDIM>::Iterator it(box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double err = std::abs((*phi)(ci) - (*w)(ci));

                err_max = std::max(err_max, err);
            }
        }
    }

    err_max = IBTK_MPI::maxReduction(err_max);

    // pout << "[PHI-W ERROR] max = " << err_max << std::endl;
}



static void
fillPhiMinusW(Pointer<PatchHierarchy<NDIM> > hierarchy, int err_idx, int phi_idx, int w_idx)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            auto err = Pointer<CellData<NDIM, double> >(patch->getPatchData(err_idx));
            auto phi = Pointer<CellData<NDIM, double> >(patch->getPatchData(phi_idx));
            auto w = Pointer<CellData<NDIM, double> >(patch->getPatchData(w_idx));

            if (err.isNull() || phi.isNull() || w.isNull()) continue;

            const Box<NDIM>& box = patch->getBox();

            for (Box<NDIM>::Iterator it(box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                (*err)(ci) = (*phi)(ci) - (*w)(ci);
            }
        }
    }
}

void
RMTMethod::postprocess()
{
    extrapolateXiNarrowBand(d_hierarchy,
                            d_phi_idx,
                            d_xi0_idx,
                            d_xi1_idx,
                            d_band_width,
                            d_num_passes);

    buildWFromXiUsingPhi0(d_hierarchy,
                          d_w_idx,
                          d_xi0_idx,
                          d_xi1_idx,
                          d_phi0_idx,
                          d_phi_idx);

    copyWToPhi(d_hierarchy, d_phi_idx, d_w_idx);

}