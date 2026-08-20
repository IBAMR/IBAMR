
#include "ibamr/RMTForceFunction.h"


#if (NDIM == 2)
#define DIV_T_FC IBAMR_FC_FUNC_(div_tmain2d, DIV_TMAIN2D)

extern "C"
{
    void DIV_T_FC(double* f11,
                  double* f12,
                  double* f21,
                  double* f22,
                  const int& f_gcw,
                  const double& GGS,
                  const double& nu_s,
                  double* Txx_s,
                  double* Txy_s,
                  double* Tyx_s,
                  double* Tyy_s,
                  const int& ts_gcw,
                  const double* phi,
                  const int& phi_gcw,
                  double* H_phi,
                  const int& H_phi_gcw,
                  const double* xi0,
                  const double* xi1,
                  const int& xi_gcw,
                  const int& ilower0,
                  const int& iupper0,
                  const int& ilower1,
                  const int& iupper1,
                  const double* dx,
                  const int& div_t_gcw,
                  double* div_t1,
                  double* div_t2);
}
#endif


#if (NDIM == 2)
void
RMTForceFunction::computePatchForce(int data_idx, SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch, double data_time)
{
    // ---------- geometry ----------
    auto geom = SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> >(patch->getPatchGeometry());
    TBOX_ASSERT(!geom.isNull());
    const double* dx = geom->getDx();
    const SAMRAI::hier::Box<NDIM>& box = patch->getBox();

    if ((d_phi_cur < 0 || d_phi_new < 0 || d_phi_scr < 0) || (d_xi0_cur < 0 || d_xi0_new < 0 || d_xi0_scr < 0) ||
        (d_xi1_cur < 0 || d_xi1_new < 0 || d_xi1_scr < 0))
    {
        pout << "[RMTForce] indices not set yet t=" << data_time << "\n";
        return;
    }

    const auto grid_geom = d_hierarchy->getGridGeometry();
    const auto& domain_box = grid_geom->getPhysicalDomain()[0];
    SAMRAI::pdat::CellIndex<NDIM> global_center((domain_box.lower() + domain_box.upper()) / 2);

    const bool do_print_patch = box.contains(global_center);
    const double eps_t = 1e-12;
    const bool do_print_time =
        do_print_patch && (!std::isfinite(d_last_print_time) || std::abs(data_time - d_last_print_time) > eps_t);

    TBOX_ASSERT(d_xi0_scr >= 0);
    TBOX_ASSERT(d_xi1_scr >= 0);
    TBOX_ASSERT(d_phi_scr >= 0);

    // ---------- output force (SideData) ----------
    auto F_data = SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> >(patch->getPatchData(data_idx));
    TBOX_ASSERT(!F_data.isNull());
    F_data->fillAll(0.0);
    // return;

    auto phi_data = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_phi_scr));
    auto xi0_data = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_xi0_scr));
    auto xi1_data = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_xi1_scr));

    // =====================================================
    // Copy CURRENT xi to SCRATCH before extrapolating xi for divT calculation
    // =====================================================
    auto phi_cur_data = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_phi_cur));
    auto xi0_cur_data = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_xi0_cur));
    auto xi1_cur_data = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_xi1_cur));

    TBOX_ASSERT(!xi0_cur_data.isNull() && !xi1_cur_data.isNull());
    TBOX_ASSERT(!xi0_data.isNull() && !xi1_data.isNull());

    // // copy CUR -> SCR (interior + ghosts that exist in SCR allocation)
    xi0_data->copy(*xi0_cur_data);
    xi1_data->copy(*xi1_cur_data);
    phi_data->copy(*phi_cur_data);

    double phi_min = 1.0e20;
    double phi_max = -1.0e20;
    int nneg = 0;

    for (SAMRAI::hier::Box<NDIM>::Iterator it(box); it; it++)
    {
        const SAMRAI::pdat::CellIndex<NDIM> ci(it());

        const double phiv = (*phi_data)(ci);

        phi_min = std::min(phi_min, phiv);
        phi_max = std::max(phi_max, phiv);

        if (phiv < 0.0) ++nneg;
    }

    // pout << "[AFTER COPY] phi_min=" << phi_min << " phi_max=" << phi_max << " nneg=" << nneg << "\n";

    //---- ensure SCRATCH ghost cells are filled (once per time) ----
    // if (!d_fill_inited)
    // {
    //     if (do_print_time) pout << "[RMTForce] ghost fill not initialized\n";
    // }
    // else
    // {
    //     const double eps = 1e-12;
    //     if (!std::isfinite(d_last_fill_time) || std::abs(data_time - d_last_fill_time) > eps)
    //     {
    //         d_last_fill_time = data_time;
    //         d_fill_op->fillData(data_time);
    //         if (do_print_time) pout << "[RMTForce] fillData DONE t=" << data_time << "\n";
    //     }
    // }

    if (do_print_time)
    {
        d_last_print_time = data_time;

        SAMRAI::pdat::CellIndex<NDIM> sidx;
        bool found_solid = false;

        auto phi_cur_data =
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_phi_cur));
        auto phi_new_data =
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_phi_new));
        auto phi_scr_data =
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_phi_scr));

        auto xi0_cur_data =
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_xi0_cur));
        auto xi0_new_data =
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_xi0_new));
        auto xi0_scr_data =
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_xi0_scr));

        auto xi1_cur_data =
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_xi1_cur));
        auto xi1_new_data =
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_xi1_new));
        auto xi1_scr_data =
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_xi1_scr));

        auto geom = Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> >(patch->getPatchGeometry());
        const double* dx = geom->getDx();
        const double margin = 1.0 * std::max(dx[0], dx[1]);

        for (SAMRAI::hier::Box<NDIM>::Iterator it(box); it; it++)
        {
            const SAMRAI::pdat::CellIndex<NDIM> ci(it());
            // if ((*phi_cur_data)(ci) < 0.0)
            if ((*phi_data)(ci) < 0.0)
            {
                sidx = ci;
                found_solid = true;
                break;
            }
        }
        // pout << "[PHI SCR DEBUG] phi_min=" << phi_min << " phi_max=" << phi_max << " nneg=" << nneg << "\n";

        if (found_solid)
        {
            auto u_data = SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> >(patch->getPatchData(data_idx));

            SAMRAI::pdat::SideIndex<NDIM> sx(sidx, 0, 0);
            SAMRAI::pdat::SideIndex<NDIM> sy(sidx, 1, 0);

            auto geom = Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> >(patch->getPatchGeometry());
            const double* dx = geom->getDx();
            const double* x_lo = geom->getXLower();
            const hier::Index<NDIM>& ilo = box.lower();

            const double x = x_lo[0] + dx[0] * (double(sidx(0) - ilo(0)) + 0.5);
            const double y = x_lo[1] + dx[1] * (double(sidx(1) - ilo(1)) + 0.5);

            const double disp_x = x - (*xi0_cur_data)(sidx);
            const double disp_y = y - (*xi1_cur_data)(sidx);
            const double disp_mag = std::sqrt(disp_x * disp_x + disp_y * disp_y);

            // pout << std::setprecision(12) << std::scientific;
            // pout << " disp x     = " << disp_x << "\n";
            // pout << " disp y     = " << disp_y << "\n";
            // pout << " disp mag   = " << disp_mag << "\n";
        }
        else
        {
            pout << "  no solid cell found on DEBUG patch\n";
        }
    }

    const int phi_gcw = phi_data->getGhostCellWidth().max();
    const int xi_gcw = xi0_data->getGhostCellWidth().max();

    const double* phi = phi_data->getPointer();
    const double* xi0 = xi0_data->getPointer();
    const double* xi1 = xi1_data->getPointer();

    // ---------- allocate work arrays used by Fortran ----------
    if (patch->getPatchData(d_F11_idx).isNull()) patch->allocatePatchData(d_F11_idx, data_time);
    if (patch->getPatchData(d_F12_idx).isNull()) patch->allocatePatchData(d_F12_idx, data_time);
    if (patch->getPatchData(d_F21_idx).isNull()) patch->allocatePatchData(d_F21_idx, data_time);
    if (patch->getPatchData(d_F22_idx).isNull()) patch->allocatePatchData(d_F22_idx, data_time);
    if (patch->getPatchData(d_Fmag_idx).isNull()) patch->allocatePatchData(d_Fmag_idx, data_time);
    if (patch->getPatchData(d_DetF_idx).isNull()) patch->allocatePatchData(d_DetF_idx, data_time);

    auto F11 = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_F11_idx));
    auto F12 = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_F12_idx));
    auto F21 = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_F21_idx));
    auto F22 = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_F22_idx));
    auto Fmag = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_Fmag_idx));
    auto DetF = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_DetF_idx));

    F11->fillAll(1.0);
    F22->fillAll(1.0);
    F12->fillAll(0.0);
    F21->fillAll(0.0);
    Fmag->fillAll(0.0);
    DetF->fillAll(1.0);

    if (patch->getPatchData(d_H_idx).isNull()) patch->allocatePatchData(d_H_idx, data_time);
    if (patch->getPatchData(d_F_idx).isNull()) patch->allocatePatchData(d_F_idx, data_time);
    if (patch->getPatchData(d_Ts_idx).isNull()) patch->allocatePatchData(d_Ts_idx, data_time);
    if (patch->getPatchData(d_div_idx).isNull()) patch->allocatePatchData(d_div_idx, data_time);

    auto H_data = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_H_idx));
    auto Ften = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_F_idx));
    auto Ts = SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(d_Ts_idx));
    auto div_data = SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> >(patch->getPatchData(d_div_idx));

    TBOX_ASSERT(!H_data.isNull() && !Ften.isNull() && !Ts.isNull() && !div_data.isNull());

    H_data->fillAll(0.0);
    Ts->fillAll(0.0);
    div_data->fillAll(0.0);

    Ften->fillAll(0.0);

    SAMRAI::hier::Box<NDIM> gbox = patch->getBox();
    gbox.grow(Ften->getGhostCellWidth());

    for (SAMRAI::hier::Box<NDIM>::Iterator it(gbox); it; it++)
    {
        const SAMRAI::pdat::CellIndex<NDIM> ci(it());

        (*Ften)(ci, 0) = 1.0;
        (*Ften)(ci, 3) = 1.0;
    }
    // ---------- pointers/GCWs (use SCRATCH phi/xi!) ----------
    double* Hphi = H_data->getPointer();
    const int H_gcw = H_data->getGhostCellWidth().max();

    double* f11 = Ften->getPointer(0);
    double* f12 = Ften->getPointer(1);
    double* f21 = Ften->getPointer(2);
    double* f22 = Ften->getPointer(3);
    const int f_gcw = Ften->getGhostCellWidth().max();

    double* Txx = Ts->getPointer(0);
    double* Txy = Ts->getPointer(1);
    double* Tyx = Ts->getPointer(2);
    double* Tyy = Ts->getPointer(3);
    const int ts_gcw = Ts->getGhostCellWidth().max();

    double* Fx_rmt = div_data->getPointer(0, 0);
    double* Fy_rmt = div_data->getPointer(1, 0);
    const int div_gcw = div_data->getGhostCellWidth().max();

    // ---------- call Fortran (compute divT into div_data) ----------
    DIV_T_FC(f11,
             f12,
             f21,
             f22,
             f_gcw,
             d_GGS,
             d_nu_s,
             Txx,
             Txy,
             Tyx,
             Tyy,
             ts_gcw,
             phi,
             phi_gcw,
             Hphi,
             H_gcw,
             xi0,
             xi1,
             xi_gcw,
             box.lower(0),
             box.upper(0),
             box.lower(1),
             box.upper(1),
             dx,
             div_gcw,
             Fx_rmt,
             Fy_rmt);

    const double margin = 1.0 * std::max(dx[0], dx[1]);

    const SAMRAI::hier::Box<NDIM>& cbox = patch->getBox();

    for (SAMRAI::hier::Box<NDIM>::Iterator it(cbox); it; it++)
    {
        const SAMRAI::pdat::CellIndex<NDIM> ci(it());

        // mask first
        if ((*phi_data)(ci) > -margin)
        {
            (*F11)(ci) = 1.0;
            (*F12)(ci) = 0.0;
            (*F21)(ci) = 0.0;
            (*F22)(ci) = 1.0;
            (*Fmag)(ci) = 0.0;
            (*DetF)(ci) = 1.0;
            continue;
        }

        const double a11 = (*Ften)(ci, 0);
        const double a12 = (*Ften)(ci, 1);
        const double a21 = (*Ften)(ci, 2);
        const double a22 = (*Ften)(ci, 3);

        (*F11)(ci) = a11;
        (*F12)(ci) = a12;
        (*F21)(ci) = a21;
        (*F22)(ci) = a22;

        const double d11 = a11 - 1.0;
        const double d22 = a22 - 1.0;

        (*Fmag)(ci) = std::sqrt(d11 * d11 + a12 * a12 + a21 * a21 + d22 * d22);
        (*DetF)(ci) = a11 * a22 - a12 * a21;
    }

    bool printed_divt_add = false;
    // ---------- add masked divT to output force ----------
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const SAMRAI::hier::Box<NDIM> sbox = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(box, axis);

        for (SAMRAI::hier::Box<NDIM>::Iterator it(sbox); it; it++)
        {
            const SAMRAI::hier::Index<NDIM>& sidx = it();
            const SAMRAI::pdat::SideIndex<NDIM> si(sidx, axis, 0);

            // Two cells adjacent to this face:
            SAMRAI::pdat::CellIndex<NDIM> c_up(sidx);
            SAMRAI::hier::IntVector<NDIM> e(0);
            e(axis) = 1;
            SAMRAI::pdat::CellIndex<NDIM> c_dn = c_up - e;

            // (Debug-safe) skip faces whose adjacent cells fall outside patch interior box
            if (!box.contains(c_up) || !box.contains(c_dn)) continue;

            const double Hup = (*H_data)(c_up); // یا (c_up,0) اگر depth>1
            const double Hdn = (*H_data)(c_dn);

            {
                const double div_val = (*div_data)(si, 0);
                const double F_before = (*F_data)(si, 0);

                {
                    (*F_data)(si, 0) += div_val;
                }

                const double F_after = (*F_data)(si, 0);

            }
        }
    }

    // ---------- gravity ramp ONLY in solid (via H_phi) ----------
    double factor = (d_Tramp > 0.0) ? (data_time / d_Tramp) : 1.0;
    factor = std::max(0.0, std::min(1.0, factor));
    const double gy = d_g0y * factor;

    const int axis = 1; // y-faces
    const SAMRAI::hier::Box<NDIM> sbox_y = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(box, axis);

    for (SAMRAI::hier::Box<NDIM>::Iterator it(sbox_y); it; it++)
    {
        const SAMRAI::hier::Index<NDIM>& sidx = it();
        SAMRAI::pdat::SideIndex<NDIM> si(sidx, axis, 0);

        SAMRAI::pdat::CellIndex<NDIM> c_up(sidx);
        SAMRAI::pdat::CellIndex<NDIM> c_dn(sidx - SAMRAI::hier::IntVector<NDIM>(0, 1));

        if (!box.contains(c_up) || !box.contains(c_dn)) continue;

        double Hface = 0.5 * ((*H_data)(c_up) + (*H_data)(c_dn));
        Hface = std::max(0.0, std::min(1.0, Hface));

        (*F_data)(si, 0) += d_rho_s * gy * Hface;
    }
}
#endif