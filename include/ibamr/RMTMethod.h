#ifndef included_ibamr_RMTMethod
#define included_ibamr_RMTMethod

#include <ibamr/config.h>
#include <ibamr/RMTForceFunction.h>

#include <ibtk/CartGridFunction.h>
#include <VisItDataWriter.h>

#include <CellVariable.h>
#include <PatchHierarchy.h>

#include <ibamr/app_namespaces.h>

#include <string>

class RMTMethod
{
public:
    RMTMethod(const std::string& object_name,
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
              double Tramp);

    Pointer<IBTK::CartGridFunction> getForceFunction() const;

    void initializeForceFunction(int phi_cur,
                                 int phi_new,
                                 int phi_scr,
                                 int xi0_cur,
                                 int xi0_new,
                                 int xi0_scr,
                                 int xi1_cur,
                                 int xi1_new,
                                 int xi1_scr);

    void initializeLevelData(double init_time);

    void registerVisItDataWriter(Pointer<VisItDataWriter<NDIM> > visit_data_writer);

    void postprocess();

private:
    std::string d_object_name;

    Pointer<PatchHierarchy<NDIM> > d_hierarchy;
    Pointer<RMTForceFunction> d_force_fcn;

    Pointer<CellVariable<NDIM, double> > d_phi0_var;
    Pointer<CellVariable<NDIM, double> > d_w_var;
    Pointer<CellVariable<NDIM, double> > d_phiw_err_var;

    int d_phi_idx = IBTK::invalid_index;
    int d_xi0_idx = IBTK::invalid_index;
    int d_xi1_idx = IBTK::invalid_index;

    int d_phi0_idx = IBTK::invalid_index;
    int d_w_idx = IBTK::invalid_index;
    int d_phiw_err_idx = IBTK::invalid_index;

    double d_band_width = 0.0;
    int d_num_passes = 0;
    int d_reinit_num_iter = 0;
    double d_reinit_dtau = 0.0;
};

#endif