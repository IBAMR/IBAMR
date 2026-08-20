#ifndef included_ibamr_RMTForceFunction
#define included_ibamr_RMTForceFunction

#include <ibamr/config.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>

#include <CellVariable.h>
#include <SideVariable.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <VariableContext.h>

#include <ibamr/app_namespaces.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>


#if (NDIM == 2)
class RMTForceFunction : public IBTK::CartGridFunction
{
public:
    void initializeGhostFill(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy)
    {
        d_hierarchy = hierarchy;

        using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;

        std::vector<ITC> comps;
        comps.emplace_back(d_phi_scr,
                           "CONSERVATIVE_LINEAR_REFINE",
                           /*use_cf_bdry_interpolation*/ false,
                           "CONSERVATIVE_COARSEN",
                           /*phys_bdry_extrap_type*/ "LINEAR",
                           /*consistent_type_2_bdry*/ false);

        comps.emplace_back(d_xi0_scr, "LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "NONE", false);

        comps.emplace_back(d_xi1_scr, "LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "NONE", false);

        d_fill_op = new IBTK::HierarchyGhostCellInterpolation();
        d_fill_op->initializeOperatorState(comps, d_hierarchy);
        d_fill_inited = true;

        pout << "[RMTForce] initializeGhostFill DONE\n";
    }

    int getF11Index() const
    {
        return d_F11_idx;
    }
    int getF12Index() const
    {
        return d_F12_idx;
    }
    int getF21Index() const
    {
        return d_F21_idx;
    }
    int getF22Index() const
    {
        return d_F22_idx;
    }
    int getFmagIndex() const
    {
        return d_Fmag_idx;
    }

    int getDetFIndex() const
    {
        return d_DetF_idx;
    }

    RMTForceFunction(const std::string& name,
                 Pointer<CellVariable<NDIM, double> > phi_var,
                 Pointer<CellVariable<NDIM, double> > xi0_var,
                 Pointer<CellVariable<NDIM, double> > xi1_var,
                 double GGS,
                 double nu_s,
                 double rho_s,
                 double g0y,
                 double Tramp)
        : d_name(name),
          d_phi_var(phi_var),
          d_xi0_var(xi0_var),
          d_xi1_var(xi1_var),
          d_GGS(GGS),
          d_nu_s(nu_s),
          d_rho_s(rho_s),
          d_g0y(g0y),
          d_Tramp(Tramp)
    {
        auto* vdb = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        d_ctx = vdb->getContext("RMT_FORCE_SCRATCH");

        const SAMRAI::hier::IntVector<NDIM> gcw_cc(1);

        d_F11_var = new SAMRAI::pdat::CellVariable<NDIM, double>(d_name + "::F11");
        d_F12_var = new SAMRAI::pdat::CellVariable<NDIM, double>(d_name + "::F12");
        d_F21_var = new SAMRAI::pdat::CellVariable<NDIM, double>(d_name + "::F21");
        d_F22_var = new SAMRAI::pdat::CellVariable<NDIM, double>(d_name + "::F22");
        d_Fmag_var = new SAMRAI::pdat::CellVariable<NDIM, double>(d_name + "::Fmag");

        d_DetF_var = new SAMRAI::pdat::CellVariable<NDIM, double>(d_name + "::DetF");

        d_F11_idx = vdb->registerVariableAndContext(d_F11_var, d_ctx, gcw_cc);
        d_F12_idx = vdb->registerVariableAndContext(d_F12_var, d_ctx, gcw_cc);
        d_F21_idx = vdb->registerVariableAndContext(d_F21_var, d_ctx, gcw_cc);
        d_F22_idx = vdb->registerVariableAndContext(d_F22_var, d_ctx, gcw_cc);
        d_Fmag_idx = vdb->registerVariableAndContext(d_Fmag_var, d_ctx, gcw_cc);

        d_DetF_idx = vdb->registerVariableAndContext(d_DetF_var, d_ctx, gcw_cc);

        d_H_var = new SAMRAI::pdat::CellVariable<NDIM, double>(d_name + "::H", 1);
        d_F_var = new SAMRAI::pdat::CellVariable<NDIM, double>(d_name + "::F", 4);
        d_Ts_var = new SAMRAI::pdat::CellVariable<NDIM, double>(d_name + "::Ts", 4);

        d_div_var = new SAMRAI::pdat::SideVariable<NDIM, double>(d_name + "::divT", 1);
        d_div_idx = vdb->registerVariableAndContext(d_div_var, d_ctx, gcw_cc);

        d_H_idx = vdb->registerVariableAndContext(d_H_var, d_ctx, gcw_cc);
        d_F_idx = vdb->registerVariableAndContext(d_F_var, d_ctx, gcw_cc);
        d_Ts_idx = vdb->registerVariableAndContext(d_Ts_var, d_ctx, gcw_cc);
    }

    void setPhi0AndWIndices(int phi0_idx, int w_idx)
    {
        d_phi0_idx = phi0_idx;
        d_w_idx = w_idx;
    }

    void setIndices(int phi_cur,
                    int phi_new,
                    int phi_scr,
                    int xi0_cur,
                    int xi0_new,
                    int xi0_scr,
                    int xi1_cur,
                    int xi1_new,
                    int xi1_scr)
    {
        d_phi_cur = phi_cur;
        d_phi_new = phi_new;
        d_phi_scr = phi_scr;
        d_xi0_cur = xi0_cur;
        d_xi0_new = xi0_new;
        d_xi0_scr = xi0_scr;
        d_xi1_cur = xi1_cur;
        d_xi1_new = xi1_new;
        d_xi1_scr = xi1_scr;
    }

    bool isTimeDependent() const override
    {
        return true;
    }

    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >) override
    {
        if (initial_time) return; // IMPORTANT: avoid calls during hierarchy init
        computePatchForce(data_idx, patch, data_time);
    }

private:
    Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_F11_var, d_F12_var, d_F21_var, d_F22_var, d_Fmag_var,
        d_DetF_var;
    int d_F11_idx = IBTK::invalid_index;
    int d_F12_idx = IBTK::invalid_index;
    int d_F21_idx = IBTK::invalid_index;
    int d_F22_idx = IBTK::invalid_index;
    int d_Fmag_idx = IBTK::invalid_index;
    int d_DetF_idx = IBTK::invalid_index;

    int d_phi0_idx = IBTK::invalid_index;
    int d_w_idx = IBTK::invalid_index;

    Pointer<PatchHierarchy<NDIM> > d_hierarchy;
    Pointer<HierarchyGhostCellInterpolation> d_fill_op;
    bool d_fill_inited = false;

    void computePatchForce(int data_idx, SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch, double data_time);

    std::string d_name;

    Pointer<CellVariable<NDIM, double> > d_phi_var, d_xi0_var, d_xi1_var;

    int d_phi_cur = -1, d_phi_new = -1, d_phi_scr = -1;
    int d_xi0_cur = -1, d_xi0_new = -1, d_xi0_scr = -1;
    int d_xi1_cur = -1, d_xi1_new = -1, d_xi1_scr = -1;

    mutable double d_last_print_time = -1.0;

    mutable double d_last_fill_time = std::numeric_limits<double>::quiet_NaN();

    double d_GGS;
    double d_nu_s;
    double d_rho_s;
    double d_g0y;
    double d_Tramp;

    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ctx;

    Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_phi_s_var, d_xi0_s_var, d_xi1_s_var;
    int d_phi_s_idx = IBTK::invalid_index;
    int d_xi0_s_idx = IBTK::invalid_index;
    int d_xi1_s_idx = IBTK::invalid_index;

    Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_H_var, d_F_var, d_Ts_var;
    int d_H_idx = IBTK::invalid_index;
    int d_F_idx = IBTK::invalid_index;
    int d_Ts_idx = IBTK::invalid_index;

    Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_div_var;
    int d_div_idx = IBTK::invalid_index;
};
#endif // NDIM == 2
#endif // included_ibamr_RMTForceFunction 