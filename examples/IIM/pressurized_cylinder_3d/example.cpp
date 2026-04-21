// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/StaggeredStokesOpenBoundaryStabilizer.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <libmesh/boundary_info.h>
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/face_tri.h>
#include <libmesh/face_tri3.h>
#include <libmesh/face_quad.h>
#include <libmesh/face_quad4.h>
#include <libmesh/gmsh_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

// Application-specific includes
//#include "FeedbackForcer.h"

namespace ModelData
{
// Elasticity model data for thin body.
// Tether (penalty) force function for thin body.
static double kappa_s_thin = 1.0;
static double eta_s_thin = 1.0;

static double L = 0.0;
static double D = 0.0;
static double H = 0.0;
static double p_e = 0.0;
static double MU = 0.0;
static double theta = 0.0;
static double zi = 0.0;
static double yi = 0.0;
static double yo = 0.0;
static double U_MAX = 0.0;
static double zo = 0.0;
static double dx = 0.0;
static double r_inner_input = 0.0;
static double r_outer_input = 0.0;
static bool use_feedback_forcer = true;
void
tether_body_force_function_thin(VectorValue<double>& F,
                                const VectorValue<double>& n,
                                const VectorValue<double>& /*N*/,
                                const TensorValue<double>& /*FF*/,
                                const libMesh::Point& x,
                                const libMesh::Point& X,
                                Elem* const elem,
                                const unsigned short /*side*/,
                                const vector<const vector<double>*>& var_data,
                                const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                                double time,
                                void* /*ctx*/)
{   
    double disp = 0.0;
    const std::vector<double>& U = *var_data[0];

    double u_bndry_n = 0.0;
    double x_kappa_n = 0.0;
    
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_kappa_n += n(d) * kappa_s_thin * (X(d) - x(d));
        u_bndry_n += n(d) * U[d];
    }

    
    for (unsigned int d = 0; d< NDIM; ++d){
        disp += (X(d) - x(d)) * (X(d) - x(d));
    }
    disp = sqrt(disp);
    double length = elem->volume();
    TBOX_ASSERT(disp < 0.75 * dx); //make sure the bdry isn't moving too much

    // The tether force is proportional to the mismatch between the positions
    // and velocities.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_s_thin * (X(d) - x(d)) + eta_s_thin * (0.0 - u_bndry_n * n(d));
    }

    //std::cout << "n is "<<n(0)<<", "<<n(1)<<", "<<n(2)<<".\n";
    //std::cout << "U is "<<U[0]<<", "<<U[1]<<", "<<U[2]<<".\n";
    //std::cout << "x current config is "<<x(0)<<", "<<x(1)<<", "<<x(2)<<".\n";
    //std::cout << "X reference config is "<<X(0)<<", "<<X(1)<<", "<<X(2)<<".\n";
    //std::cout << "F tether force is "<<F(0)<<", "<<F(1)<<", "<<F(2)<<".\n";
    //if(time != 2e-6){
    //    std::cin.get();
    //}
    return;
}

inline unsigned int
idx(const unsigned int nr, const unsigned int i, const unsigned int j)
{
    return i + j * nr;

    return libMesh::invalid_uint;
}

static double x_loc_min, x_loc_max, z_loc_min, z_loc_max, y_loc_min, y_loc_max;
} // namespace ModelData
using namespace ModelData;

inline double
cos_shape_fcn(const double x, const double radius)
{
    if (std::abs(x) > radius)
        return 0.0;
    else
        return 0.5 * (1.0 + cos(M_PI * x / radius));
} // cos_shape_fcn
inline double
plug_flow_profile(const IBTK::Point& x,
               const IBTK::Point& x_center,
               const double R,
               const double epsilon,
               const int bdry_normal_axis)
{
    double r_sq = 0.0;
    for (int d = 0; d < NDIM; ++d) r_sq += (d == bdry_normal_axis) ? 0.0 : (x[d] - x_center[d]) * (x[d] - x_center[d]);
    const double r = sqrt(r_sq);
    if (r > R)
        return 0.0;
    else if (r < R - epsilon)
        return 1.0;
    else
        return cos_shape_fcn(r - (R - epsilon), epsilon);
} // plug_flow_profile

inline double
parabolic_flow_profile(const IBTK::Point& x,
                       const IBTK::Point& x_center,
                       const double R,
                       const int bdry_normal_axis)
{
    double r_sq = 0.0;
    for (int d = 0; d < NDIM; ++d)
        r_sq += (d == bdry_normal_axis) ? 0.0 : (x[d] - x_center[d]) * (x[d] - x_center[d]);
    const double r = sqrt(r_sq);

    if (r > R)
        return 0.0;  // Outside the pipe radius
    else
        return 1.0 - (r_sq / (R * R));  // Parabolic profile normalized to peak at 1.0
}


inline double
normal_profile(const double x, const double x_center, const double H, const double epsilon)
{
    const double h = abs(x - x_center);
    if (h > H)
        return 0.0;
    else if (h < H - epsilon)
        return 1.0;
    else
        return cos_shape_fcn(h - (H - epsilon), epsilon);
} // normal_profile

// Bessel function of the first kind of order 0
std::complex<double> besselJ0(std::complex<double> z) {
    std::complex<double> sum(1.0, 0.0);
    std::complex<double> term(1.0, 0.0);
    for (int k = 1; k <= 50; ++k) 
    {
        term *= -0.25 * z * z / static_cast<double>(k * k);
        sum += term;
        if (std::abs(term) < 1e-10 * std::abs(sum)) break;
    }
    return sum;
}

double womersley_profile(double R, double omega, double phi0, double rho, double nu, double dpdz, double r, double t, double epsilon)
{
    if (r > R) 
        return 0.0;


    double alpha = R * sqrt(omega / nu);  // Womersley number
    std::complex<double> j(0.0, 1.0); 

    std::complex<double> lambda_c = alpha * std::pow(j, 3.0/2.0);
    
    std::complex<double> A = dpdz / (j * rho * omega);

    std::complex<double> t1 = besselJ0(lambda_c*r/R);
    std::complex<double> t2 = besselJ0(lambda_c);
    std::complex<double> t3 =  std::exp(j * (omega * t - phi0));
    
    std::complex<double> velocity = A * (1 - t1 / (t2 + 1e-30)) * t3;

    return velocity.real();

}


inline Box<NDIM>
get_domain_bdry_box(const double H,
                    const int axis,
                    const int side,
                    Pointer<CartesianGridGeometry<NDIM> > grid_geom,
                    const IntVector<NDIM>& ratio)
{
    const auto& domain_box = Box<NDIM>::refine(grid_geom->getPhysicalDomain()[0], ratio);
    std::array<double, NDIM> dx;
    for (int d = 0; d < NDIM; ++d) dx[d] = grid_geom->getDx()[d] / static_cast<double>(ratio(d));
    Box<NDIM> bdry_box = domain_box;
    const int offset = std::ceil(H / dx[axis]);
    const bool at_lower_bdry = side == 0;
    if (at_lower_bdry)
    {
        bdry_box.upper(axis) = domain_box.lower(axis) + offset;
    }
    else
    {
        bdry_box.lower(axis) = domain_box.upper(axis) - offset;
    }
    return bdry_box;
} // get_domain_bdry_box

inline Box<NDIM>
get_outlet_bdry_box(const CellIndex<NDIM>& i_center,
                    const double R,
                    const double H,
                    const int axis,
                    const int side,
                    Pointer<CartesianGridGeometry<NDIM> > grid_geom,
                    const IntVector<NDIM>& ratio)
{
    std::array<double, NDIM> dx;
    for (int d = 0; d < NDIM; ++d) dx[d] = grid_geom->getDx()[d] / static_cast<double>(ratio(d));
    Box<NDIM> outlet_box(i_center, i_center);
    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis)
        {
            const int offset = std::ceil(H / dx[axis]);
            const bool at_lower_bdry = side == 0;
            if (at_lower_bdry)
                outlet_box.upper(d) += offset;
            else
                outlet_box.lower(d) -= offset;
        }
        else
        {
            const int offset = std::ceil(R / dx[axis]);
            outlet_box.upper(d) += offset;
            outlet_box.lower(d) -= offset;
        }
    }
    return outlet_box * get_domain_bdry_box(H, axis, side, grid_geom, ratio);
} // get_outlet_bdry_box

 // namespace

 class ParabolicFeedbackForcer : public CartGridFunction
 {
 public:
     ParabolicFeedbackForcer(const INSHierarchyIntegrator* fluid_solver, const Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
         : d_fluid_solver(fluid_solver), d_patch_hierarchy(patch_hierarchy)
     {
     } // ParabolicFeedbackForcer
 
     ~ParabolicFeedbackForcer() = default;
 
     bool isTimeDependent() const
     {
         return true;
     } // isTimeDependent
 
     void setDataOnPatchHierarchy(int data_idx,
                                  Pointer<hier::Variable<NDIM> > /*var*/,
                                  Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                  double /*data_time*/,
                                  bool initial_time,
                                  int coarsest_ln_in,
                                  int finest_ln_in)
     {
         // We wish to apply constraint forces of the form f = kappa * (U0 * u_profile - u_grid) that force the flow to
         // assume a particular profile near inlets and outlets. To avoid adding or removing momentum, we want the net
         // force associated with the constraint force to be zero. (Perhaps we should also modify it to have zero net
         // torque.) We also do not know the actual flow velocity near the boundary and do not want to force it to take a
         // particular value --- only to adopt a particular overall profile.
         //
         // The function proceeds in two steps:
         //
         // 1) Determine U0 so that constraint force has net zero.
         //
         // 2) Apply the force.
         //
         // We cannot perform the two steps at the same time because competing step one could involve contributions from
         // multiple processors.
 
         const int coarsest_hier_ln = 0;
         const int finest_hier_ln = patch_hierarchy->getFinestLevelNumber();
         const int coarsest_ln = (coarsest_ln_in == invalid_level_number ? coarsest_hier_ln : coarsest_ln_in);
         const int finest_ln = (finest_ln_in == invalid_level_number ? finest_hier_ln : finest_ln_in);
         Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_hierarchy->getGridGeometry();
 
         // Setup regions where forces are applied.
         const double* const dx_coarsest = grid_geom->getDx();
         const IntVector<NDIM>& finest_ratio = patch_hierarchy->getPatchLevel(finest_hier_ln)->getRatio();
         std::array<double, NDIM> dx_finest;
         for (int d = 0; d < NDIM; ++d) dx_finest[d] = dx_coarsest[d] / static_cast<double>(finest_ratio(d));
         const double epsilon = 1.0 * dx_finest[0];
         for (int d = 0; d < NDIM; ++d) TBOX_ASSERT(epsilon == 1.0 * dx_finest[d]);
         const double H = 2.0 * dx_finest[0];
 
         const int side_wgt_idx = d_fluid_solver->getHierarchyMathOps()->getSideWeightPatchDescriptorIndex();
 
         // Position information for the outlets.
         const int n_outlets = 2;
         const std::vector<IBTK::Point> x_outlet = { IBTK::Point(0.0, 2.5, 2.5), IBTK::Point(5.0, 2.5, 2.5)  };
         const std::vector<int> outlet_axis = { 0, 0 };
         const std::vector<int> outlet_side = { 0, 1 };
	 
	 std::cout << "R inner is "<< r_inner_input<<", r outer is "<< r_outer_input<<"\n\n"; 
         const std::vector<double> r_outlet_inner = { r_inner_input, r_inner_input};
         const std::vector<double> r_outlet_outer = { r_outer_input, r_outer_input};
         //const std::vector<double> r_outlet_inner = { 1.0 - 0 * dx_finest[0], 1.0 - 0 * dx_finest[0]};
         //const std::vector<double> r_outlet_outer = { 1.0 + 0.0 * dx_finest[0], 1.0 + 0.0 * dx_finest[0]};
         // Setup the penalty parameter.
         const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
         const double dt = d_fluid_solver->getCurrentTimeStepSize();
         const double rho = d_fluid_solver->getStokesSpecifications()->getRho();
         const double kappa = cycle_num >= 0 ? 0.125 * rho / dt : 0.0;
 
         // Compute scale factors for velocity profiles to ensure net forces are zero.
         std::vector<IBTK::Vector> u_profile_integral(n_outlets, IBTK::Vector::Zero());
         std::vector<IBTK::Vector> u_grid_integral(n_outlets, IBTK::Vector::Zero());
         for (int ln = coarsest_hier_ln; ln <= finest_hier_ln; ++ln)
         {
             Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
             const IntVector<NDIM>& ratio = level->getRatio();
             for (PatchLevel<NDIM>::Iterator p(level); p; p++)
             {
                 Pointer<Patch<NDIM> > patch = level->getPatch(p());
                 const Box<NDIM>& patch_box = patch->getBox();
 
                 Pointer<SideData<NDIM, double> > U_current_data =
                     patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
                 Pointer<SideData<NDIM, double> > U_new_data =
                     patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());
                 TBOX_ASSERT(U_current_data);
 
                 Pointer<SideData<NDIM, double> > side_wgt_data = patch->getPatchData(side_wgt_idx);
 
                 for (int k = 0; k < n_outlets; ++k)
                 {
                     const double R_inner = r_outlet_inner[k];
                     const double R_outer = r_outlet_outer[k];
                     const int axis = outlet_axis[k];
                     const int side = outlet_side[k];
                     const bool at_lower_bdry = side == 0;
                     const double x_bdry = at_lower_bdry ? grid_geom->getXLower()[axis] : grid_geom->getXUpper()[axis];
                     const IBTK::Point& x_center = x_outlet[k];
                     const CellIndex<NDIM> i_center = IndexUtilities::getCellIndex(x_center, grid_geom, ratio);
 
                     // Determine the region where the forcing should be applied.
                     const Box<NDIM> bdry_box = get_outlet_bdry_box(i_center, R_inner, H, axis, side, grid_geom, ratio);
                     const Box<NDIM> patch_bdry_box = bdry_box * patch_box;
                     if (patch_bdry_box.empty()) continue;
 
                     // Sum up the local contributions.
                     for (int d = 0; d < NDIM; ++d)
                     {
                         for (SideIterator<NDIM> it(patch_bdry_box, d); it; it++)
                         {
                             const SideIndex<NDIM>& i_s = it();
                             const auto x = IndexUtilities::getSideCenter(*patch, i_s);
                             double r_sq = 0.0;
                             for (int d = 0; d < NDIM; ++d)
                                 r_sq += (d == axis) ? 0.0 : std::pow(x(d) - x_center(d), 2.0);
                             const double r = std::sqrt(r_sq);
 
                             // We can skip locations that do not involve staggered DOFs.
                             const double dx = (*side_wgt_data)(i_s);
                             if (dx == 0.0) continue;
 
                             // Accumulate local contributions.
                             const double wgt = normal_profile(x(axis), x_bdry, H, epsilon) *
                                                parabolic_flow_profile(x, x_center, R_outer, axis);
                             if (r <= R_inner && wgt > 0.0 && R_inner != 0.0)
                             {
                                 const double u_profile = parabolic_flow_profile(x, x_center, R_outer, axis);
                                 u_profile_integral[k](d) += u_profile * wgt * dx;
 
                                 const double u_grid_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                                 const double u_grid_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                                 const double u_grid =
                                     (cycle_num > 0) ? 0.5 * (u_grid_new + u_grid_current) : u_grid_current;
                                 u_grid_integral[k](d) += u_grid * wgt * dx;
                             }
                         }
                     }
                 }
             }
         }
 
         // Find the profile that will apply zero total force.
         //
         // We want: int (U0 * u_profile - u_grid) wgt dx = 0 ===> U0 = int u_grid wgt dx / int u_profile wgt dx.
         std::vector<IBTK::Vector> U0_profile(n_outlets, IBTK::Vector::Zero());
         for (int k = 0; k < n_outlets; ++k)
         {
             SAMRAI_MPI::sumReduction(u_profile_integral[k].data(), NDIM);
             SAMRAI_MPI::sumReduction(u_grid_integral[k].data(), NDIM);
             for (int d = 0; d < NDIM; ++d) {
                U0_profile[k](d) = u_grid_integral[k](d) / u_profile_integral[k](d);
                //std::cout << "k = "<< k<< "d = "<< d <<".\n";
                //std::cout<<"a * u0 = "<<U0_profile[k](d) * parabolic_flow_profile(x, x_center, R_outer, axis)<<".\n\n";
            }

         }
 
         // Apply the force to the grid on the specified levels of the patch hierarchy (which could be a subset of the
         // full range of levels in the patch hierarchy).
         std::vector<IBTK::Vector> u_constraint_net_force(n_outlets, IBTK::Vector::Zero());
         for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
         {
             Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
             const IntVector<NDIM>& ratio = level->getRatio();
             const auto& domain_box = Box<NDIM>::refine(grid_geom->getPhysicalDomain()[0], ratio);
             for (PatchLevel<NDIM>::Iterator p(level); p; p++)
             {
                 Pointer<Patch<NDIM> > patch = level->getPatch(p());
                 const Box<NDIM>& patch_box = patch->getBox();
 
                 Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(data_idx);
                 TBOX_ASSERT(F_data);
                 F_data->fillAll(0.0);
                 if (initial_time) continue;
 
                 Pointer<SideData<NDIM, double> > U_current_data =
                     patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
                 Pointer<SideData<NDIM, double> > U_new_data =
                     patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());
                 TBOX_ASSERT(U_current_data);
 
                 Pointer<SideData<NDIM, double> > side_wgt_data = patch->getPatchData(side_wgt_idx);
 
                 for (int k = 0; k < n_outlets; ++k)
                 {
                     const double R_inner = r_outlet_inner[k];
                     const double R_outer = r_outlet_outer[k];
                     const int axis = outlet_axis[k];
                     const int side = outlet_side[k];
                     const bool at_lower_bdry = side == 0;
                     const double x_bdry = at_lower_bdry ? grid_geom->getXLower()[axis] : grid_geom->getXUpper()[axis];
                     const IBTK::Point& x_center = x_outlet[k];
                     const CellIndex<NDIM> i_center = IndexUtilities::getCellIndex(x_center, grid_geom, ratio);
 
                     // Determine the region where the forcing should be applied.
                     const Box<NDIM> bdry_box = get_outlet_bdry_box(i_center, R_outer, H, axis, side, grid_geom, ratio);
                     const Box<NDIM> patch_bdry_box = bdry_box * patch_box;
                     if (patch_bdry_box.empty()) continue;
 
                     // Apply the local contributions.
                     for (int d = 0; d < NDIM; ++d)
                     {
                         for (SideIterator<NDIM> it(patch_bdry_box, d); it; it++)
                         {
                             const SideIndex<NDIM>& i_s = it();
                             const auto x = IndexUtilities::getSideCenter(*patch, i_s);
                             double r_sq = 0.0;
                             for (int d = 0; d < NDIM; ++d)
                                 r_sq += (d == axis) ? 0.0 : std::pow(x(d) - x_center(d), 2.0);
                             const double r = std::sqrt(r_sq);
 
                             // We can skip locations that do not involve staggered DOFs.
                             const double dx = (*side_wgt_data)(i_s);
                             if (dx == 0.0) continue;
 
                             // Apply the local contributions.
                             const double wgt = normal_profile(x(axis), x_bdry, H, epsilon) *
                                                parabolic_flow_profile(x, x_center, R_outer, axis);
                             if (r <= R_outer && wgt > 0.0 && R_outer != 0.0)
                             {
                                 const double u_grid_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                                 const double u_grid_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                                 const double u_grid =
                                     (cycle_num > 0) ? 0.5 * (u_grid_new + u_grid_current) : u_grid_current;
                                 const double u_profile =
                                     U0_profile[k](d) * parabolic_flow_profile(x, x_center, R_outer, axis);
                                 if (r <= R_inner)
                                 {
                                     const double F = kappa * (u_profile - u_grid) * wgt;
                                     u_constraint_net_force[k](d) += F * dx;
                                     (*F_data)(i_s) += F;
                                 }
                             }
                         }
                     }
                 }
             }
         }
 
         // Confirm the the force near each outlet has mean zero.
         for (int k = 0; k < n_outlets; ++k)
         {
             SAMRAI_MPI::sumReduction(u_constraint_net_force[k].data(), NDIM);
             plog << "net constraint force for boundary " << k << ": (" << u_constraint_net_force[k](0) << ","
                  << u_constraint_net_force[k](1) << ")\n";
             for (int d = 0; d < NDIM; ++d)
             {
                 if (!IBTK::abs_equal_eps(u_constraint_net_force[k](d), 0.0))
                 {
                     pout << "WARNING: at boundary " << k << ", component " << d
                          << " of the constraint force has nonzero mean!\n";
                 }
             }
         }
     }
 
     void setDataOnPatch(const int data_idx,
                         Pointer<hier::Variable<NDIM> > /*var*/,
                         Pointer<Patch<NDIM> > patch,
                         const double /*data_time*/,
                         const bool initial_time,
                         Pointer<PatchLevel<NDIM> > /*patch_level*/)
     {
         Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(data_idx);
         TBOX_ASSERT(F_data);
         F_data->fillAll(0.0);
         if (initial_time) return;
         TBOX_ERROR("implementation requires collective operations: use setDataOnPatchHierarchy instead\n");
         return;
     } // setDataOnPatch
 
 private:
     ParabolicFeedbackForcer() = delete;
 
     ParabolicFeedbackForcer(const ParabolicFeedbackForcer& from) = delete;
 
     ParabolicFeedbackForcer& operator=(const ParabolicFeedbackForcer& that) = delete;
 
     const INSHierarchyIntegrator* const d_fluid_solver;
     Pointer<PatchHierarchy<NDIM> > d_patch_hierarchy;
 };


static ofstream max_displacement_stream,flow_rate_stream;
// Function prototypes
void compute_velocity_profile(tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                              const int u_idx,
                              const double /*data_time*/,
                              const string& /*data_dump_dirname*/);

void compute_pressure_profile(tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                              const int p_idx,
                              const double /*data_time*/,
                              const string& /*data_dump_dirname*/);

void postprocess_data(tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      tbox::Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
                      const int iteration_num,
                      const double loop_time,
                      const string& data_dump_dirname);
void postprocess_displacement_data(MeshBase &mesh, System &dX_system);

void output_data(tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 tbox::Pointer<IBHierarchyIntegrator> ins_integrator,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

void compute_flow_rate(const tbox::Pointer<PatchHierarchy<NDIM> > hierarchy,
                       const int U_idx,
                       const double loop_time,
                       const int wgt_sc_idx);
/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-8");
    PetscOptionsSetValue(nullptr, "-stokes_ksp_atol", "1e-8");

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        tbox::Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        tbox::Pointer<tbox::Database> input_db = app_initializer->getInputDatabase();
        r_inner_input = input_db->getDouble("RADIUS_INNER");
        r_outer_input = input_db->getDouble("RADIUS_OUTER");
        use_feedback_forcer = input_db->getBool("USE_FEEDBACK_FORCER");
	// Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string thin_exodus_filename = app_initializer->getExodusIIFilename("cylinder_thin");

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();
        const bool USE_TANGENTIAL_TREATMENT =  input_db->getBoolWithDefault("USE_TANGENTIAL_TREATMENT", false);
        // Create FE mesh.
        D = input_db->getDouble("D");              // channel height (cm)
        L = input_db->getDouble("L");              // channel length (cm)
        H = input_db->getDouble("H");              // channel length (cm)
        const double H = input_db->getDouble("H"); // wall thickness (cm)

        zo = input_db->getDouble("Zo");
        yo = input_db->getDouble("Yo");
        yi = input_db->getDouble("Yi");
        zi = input_db->getDouble("Zi");
        U_MAX = input_db->getDouble("U_MAX");
         dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        MU = input_db->getDouble("MU");
        theta = input_db->getDouble("THETA");
        p_e = input_db->getDouble("P_E");

        Mesh cylinder_mesh_thin(init.comm(), NDIM - 1);
        string cyl_mesh_filename = input_db->getString("CYL_MESH_FILENAME");

        cylinder_mesh_thin.get_boundary_info().clear_boundary_node_ids();
        libMesh::GmshIO gmsh_io(cylinder_mesh_thin);
        gmsh_io.read(cyl_mesh_filename);
        BoundaryInfo& boundary_info = cylinder_mesh_thin.get_boundary_info();
        /*
        const unsigned int NXi_elem = ceil(L / ds);
        const unsigned int NRi_elem = ceil(M_PI * D / ds);
        int node_id = 0;
        cylinder_mesh_thin.reserve_nodes(NRi_elem * (NXi_elem + 1));
        cylinder_mesh_thin.reserve_elem(NRi_elem * NXi_elem);

        for (unsigned int j = 0; j <= NXi_elem; j++)
        {
            for (unsigned int i = 0; i <= NRi_elem - 1; i++)
            {
                const double theta = 2.0 * M_PI * static_cast<Real>(i) / static_cast<Real>(NRi_elem);
                cylinder_mesh_thin.add_point(libMesh::Point(L * static_cast<Real>(j) / static_cast<Real>(NXi_elem),
                                                            0.5 * H + 0.5 * D * cos(theta),
                                                            0.5 * H + 0.5 * D * sin(theta)),
                                             node_id++);
            }
        }

        for (unsigned int j = 0; j <= NXi_elem - 1; j++)
        {
            for (unsigned int i = 0; i <= NRi_elem - 2; i++)
            {
                Elem* elem = cylinder_mesh_thin.add_elem(new Quad4);
                elem->set_node(0) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, i, j));
                elem->set_node(1) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, i + 1, j));
                elem->set_node(2) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, i + 1, j + 1));
                elem->set_node(3) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, i, j + 1));
            }
        }

        for (unsigned int j = 0; j <= NXi_elem - 1; j++)
        {
            Elem* elem = cylinder_mesh_thin.add_elem(new Quad4);
            elem->set_node(0) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, NRi_elem - 1, j));
            elem->set_node(1) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, 0, j));
            elem->set_node(2) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, 0, j + 1));
            elem->set_node(3) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, NRi_elem - 1, j + 1));
        }

        MeshBase::const_element_iterator el_end = cylinder_mesh_thin.elements_end();
        for (MeshBase::const_element_iterator el = cylinder_mesh_thin.elements_begin(); el != el_end; ++el)
        {
            const Elem* elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (at_mesh_bdry)
                {
                    BoundaryInfo& boundary_info_cylinder = cylinder_mesh_thin.get_boundary_info();
                    boundary_info_cylinder.add_side(elem, side, FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID);
                }
            }
        }
        */
        cylinder_mesh_thin.prepare_for_use();

        vector<MeshBase*> meshes(1);
        meshes[0] = &cylinder_mesh_thin;

        // Setup data for imposing constraints.
        kappa_s_thin = input_db->getDoubleWithDefault("KAPPA_S_THIN", kappa_s_thin);
        eta_s_thin = input_db->getDoubleWithDefault("ETA_S_THIN", eta_s_thin);

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        tbox::Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        tbox::Pointer<IIMethod> ib_method_ops =
            new IIMethod("IIMethod",
                         app_initializer->getComponentDatabase("IIMethod"),
                         meshes,
                         app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        tbox::Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        tbox::Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy =
            new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        tbox::Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        tbox::Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        tbox::Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        tbox::Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IBFE solver.
        ib_method_ops->initializeFEEquationSystems();
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> sys_data(1);
        sys_data[0] = SystemData(IIMethod::VELOCITY_SYSTEM_NAME, vars);
        IIMethod::LagSurfaceForceFcnData tether_body_force_thin_data(tether_body_force_function_thin, sys_data);

        //~ ib_method_ops->registerLagBodyForceFunction(tether_body_force_data, 0);
        ib_method_ops->registerLagSurfaceForceFunction(tether_body_force_thin_data, 0);
        if (USE_TANGENTIAL_TREATMENT)
           ib_method_ops->registerTangentialVelocityMotion(0);
        ib_method_ops->initializeFEEquationSystems();

        // Create Eulerian initial condition specification objects.

        tbox::Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);

        tbox::Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

        // Create Eulerian boundary condition specification objects.
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL));
        const bool periodic_domain = grid_geometry->getPeriodicShift().min() > 0;
        if (!periodic_domain)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();
                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();
                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }

            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

            if (input_db->keyExists("BoundaryStabilization"))
            {
                time_integrator->registerBodyForceFunction(new StaggeredStokesOpenBoundaryStabilizer(
                    "BoundaryStabilization",
                    app_initializer->getComponentDatabase("BoundaryStabilization"),
                    navier_stokes_integrator,
                    grid_geometry));
            }
        }
        // Create Eulerian body force function specification objects.

        //time_integrator->registerBodyForceFunction(new FeedbackForcer(navier_stokes_integrator, patch_hierarchy));
	if(use_feedback_forcer){
            time_integrator->registerBodyForceFunction(new ParabolicFeedbackForcer(navier_stokes_integrator, grid_geometry));
	}

        // Set up visualization plot file writers.
        tbox::Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_io_thin(uses_exodus ? new ExodusII_IO(cylinder_mesh_thin) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Set up locations to get velocity data.
        x_loc_min = input_db->getDoubleWithDefault("X_LOC_MIN", 0.1 * L);
        x_loc_max = input_db->getDoubleWithDefault("X_LOC_MAX", 0.9 * L);
        z_loc_min = input_db->getDoubleWithDefault("Z_LOC_MIN", 0.0);
        y_loc_min = input_db->getDoubleWithDefault("Y_LOC_MIN", 0.0);
        y_loc_max = input_db->getDoubleWithDefault("Y_LOC_MAX", L);
        z_loc_max = input_db->getDoubleWithDefault("Z_LOC_MAX", L);

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        EquationSystems* equation_systems_thin = ib_method_ops->getFEDataManager(0)->getEquationSystems();
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                exodus_io_thin->write_timestep(
                    thin_exodus_filename, *equation_systems_thin, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        if (SAMRAI_MPI::getRank() == 0)
        {
            max_displacement_stream.open("max_displacement");
            flow_rate_stream.open("flow_rate_outflow_kappa_" + std::to_string(int(kappa_s_thin)) + ".curve",
                                  ios_base::out | ios_base::trunc);
            flow_rate_stream.precision(10);
        }
        // Main time step loop.

        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";
            tbox::Pointer<hier::Variable<NDIM> > u_var = time_integrator->getVelocityVariable();
            tbox::Pointer<VariableContext> current_ctx = time_integrator->getCurrentContext();

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                    exodus_io_thin->write_timestep(
                        thin_exodus_filename, *equation_systems_thin, iteration_num / viz_dump_interval + 1, loop_time);

                    exodus_io_thin->write_timestep(
                        thin_exodus_filename, *equation_systems_thin, iteration_num / viz_dump_interval + 1, loop_time);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }

            {
                VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
                tbox::Pointer<hier::Variable<NDIM> > u_var = time_integrator->getVelocityVariable();
                const tbox::Pointer<VariableContext> u_ctx = time_integrator->getCurrentContext();
                const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
                const int coarsest_ln = 0;
                const int finest_ln = patch_hierarchy->getFinestLevelNumber();
                HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
                hier_math_ops.setPatchHierarchy(patch_hierarchy);
                hier_math_ops.resetLevels(coarsest_ln, finest_ln);
                const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

                compute_flow_rate(patch_hierarchy, u_idx, loop_time, wgt_sc_idx);
            }
            pout << "\nWriting state data...\n\n";

            postprocess_data(patch_hierarchy,
                             navier_stokes_integrator,
                             cylinder_mesh_thin,
                             equation_systems_thin,
                             iteration_num,
                             loop_time,
                             postproc_data_dump_dirname);
            postprocess_displacement_data(equation_systems_thin->get_mesh(),equation_systems_thin->get_system(IIMethod::COORD_MAPPING_SYSTEM_NAME));
        }

        // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        tbox::Pointer<hier::Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
        const tbox::Pointer<VariableContext> u_ctx = navier_stokes_integrator->getCurrentContext();
        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
        const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();

        const tbox::Pointer<hier::Variable<NDIM> > p_var = time_integrator->getPressureVariable();
        const tbox::Pointer<VariableContext> p_ctx = time_integrator->getCurrentContext();

        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);
        const int p_cloned_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_cloned_idx, loop_time);
        }
        u_init->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);
        p_init->setDataOnPatchHierarchy(p_cloned_idx, p_var, patch_hierarchy, loop_time - 0.5 * dt);

        compute_pressure_profile(patch_hierarchy, p_idx, loop_time, postproc_data_dump_dirname);
        compute_velocity_profile(patch_hierarchy, u_idx, loop_time, postproc_data_dump_dirname);

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_sc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
        pout << std::setprecision(16) << "Error in u at time " << loop_time << ":\n"
             << "  L1-norm:  " << hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx) << "\n"
             << "  L2-norm:  " << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << "\n"
             << "  max-norm: " << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.subtract(p_idx, p_idx, p_cloned_idx);
        pout << "Error in p at time " << loop_time - 0.5 * dt << ":\n"
             << "  L1-norm:  " << hier_cc_data_ops.L1Norm(p_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(p_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(p_idx, wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        if (dump_viz_data && uses_visit)
        {
            navier_stokes_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);

            exodus_io_thin->write_timestep(
                thin_exodus_filename, *equation_systems_thin, iteration_num / viz_dump_interval + 1, loop_time);
        }

        if (SAMRAI_MPI::getRank() == 0)
        {
            max_displacement_stream.close();
            flow_rate_stream.close();
            
        }
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
} // main

void
compute_velocity_profile(tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                         const int u_idx,
                         const double /*data_time*/,
                         const string& /*data_dump_dirname*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
    hier_math_ops.resetLevels(finest_ln, finest_ln);
    const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
    const double X_min[3] = { 0.1 * L, y_loc_min, z_loc_min };
    const double X_max[3] = { 0.9 * L, y_loc_max, y_loc_max };
    int qp_tot = 0;
    double u_Eulerian_L2_norm = 0.0;
    double u_Eulerian_max_norm = 0.0;
    vector<double> pos_values;
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        tbox::Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            tbox::Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();

            const tbox::Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();

            const double* const patch_dx = patch_geom->getDx();

            // Entire box containing the required data.
            Box<NDIM> box(IndexUtilities::getCellIndex(
                              &X_min[0], patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper),
                          IndexUtilities::getCellIndex(
                              &X_max[0], patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper));
            // Part of the box on this patch
            Box<NDIM> trim_box = patch_box * box;
            BoxList<NDIM> iterate_box_list = trim_box;

            // Trim the box covered by the finer region
            BoxList<NDIM> covered_boxes;
            if (ln < finest_ln)
            {
                BoxArray<NDIM> refined_region_boxes;
                tbox::Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
                refined_region_boxes = next_finer_level->getBoxes();
                refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM> refined_box = refined_region_boxes[i];
                    const Box<NDIM> covered_box = trim_box * refined_box;
                    covered_boxes.unionBoxes(covered_box);
                }
            }
            iterate_box_list.removeIntersections(covered_boxes);

            // Loop over the boxes and store the location and interpolated value.
            tbox::Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
            const tbox::Pointer<CellData<NDIM, double> > wgt_cc_data = patch->getPatchData(wgt_cc_idx);

            for (BoxList<NDIM>::Iterator lit(iterate_box_list); lit; lit++)
            {
                const Box<NDIM>& iterate_box = *lit;
                for (Box<NDIM>::Iterator bit(iterate_box); bit; bit++)
                {
                    const CellIndex<NDIM>& lower_idx = *bit;

                    const double yu = patch_x_lower[1] + patch_dx[1] * (lower_idx(1) - patch_lower(1) + 0.5);
                    const double zu = patch_x_lower[2] + patch_dx[2] * (lower_idx(2) - patch_lower(2) + 0.5);
                    const double xu = patch_x_lower[0] + patch_dx[0] * (lower_idx(0) - patch_lower(0));

                    double u_ex, v_ex, w_ex;

                    if (sqrt((yu - 0.5 * H) * (yu - 0.5 * H) + (zu - 0.5 * H) * (zu - 0.5 * H)) >= D / 2)
                    {
                        u_ex = 0.0;
                        v_ex = 0.0;
                        w_ex = 0.0;
                    }
                    else
                    {
                        u_ex =
                            U_MAX *
                            (1.0 - 4.0 * ((yu - 0.5 * H) * (yu - 0.5 * H) + (zu - 0.5 * H) * (zu - 0.5 * H)) / (D * D));
                        v_ex = 0.0;
                        w_ex = 0.0;
                    }

                    const double u0 = (*u_data)(SideIndex<NDIM>(lower_idx, 0, SideIndex<NDIM>::Lower));
                    const double v0 = (*u_data)(SideIndex<NDIM>(lower_idx, 1, SideIndex<NDIM>::Lower));
                    const double w0 = (*u_data)(SideIndex<NDIM>(lower_idx, 2, SideIndex<NDIM>::Lower));

                    if (xu > 0.4 * L && xu < 0.6 * L)
                    {
                        qp_tot += 1;
                        u_Eulerian_L2_norm += std::abs(u0 - u_ex) * std::abs(u0 - u_ex) * (*wgt_cc_data)(lower_idx);
                        u_Eulerian_L2_norm += std::abs(v0 - v_ex) * std::abs(v0 - v_ex) * (*wgt_cc_data)(lower_idx);
                        u_Eulerian_L2_norm += std::abs(w0 - w_ex) * std::abs(w0 - w_ex) * (*wgt_cc_data)(lower_idx);

                        u_Eulerian_max_norm = std::max(u_Eulerian_max_norm, std::abs(u0 - u_ex));
                        u_Eulerian_max_norm = std::max(u_Eulerian_max_norm, std::abs(v0 - v_ex));
                        u_Eulerian_max_norm = std::max(u_Eulerian_max_norm, std::abs(w0 - w_ex));
                    }
                }
            }
        }
    }

    SAMRAI_MPI::sumReduction(&qp_tot, 1);
    SAMRAI_MPI::sumReduction(&u_Eulerian_L2_norm, 1);
    SAMRAI_MPI::maxReduction(&u_Eulerian_max_norm, 1);
    u_Eulerian_L2_norm = sqrt(u_Eulerian_L2_norm);

    pout << " u_Eulerian_L2_norm = " << u_Eulerian_L2_norm << "\n\n";
    pout << " u_Eulerian_max_norm = " << u_Eulerian_max_norm << "\n\n";

    return;
} // compute_velocity_profile

void
output_data(tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            tbox::Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    tbox::Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();
    return;
} // output_data

void
postprocess_data(tbox::Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                 tbox::Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                 Mesh& mesh,
                 EquationSystems* equation_systems,
                 const int /*iteration_num*/,
                 const double loop_time,
                 const string& /*data_dump_dirname*/)
{
    const unsigned int dim = mesh.mesh_dimension();
    double F_integral[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;

    System& x_system = equation_systems->get_system(IIMethod::COORDS_SYSTEM_NAME);
    System& U_system = equation_systems->get_system(IIMethod::VELOCITY_SYSTEM_NAME);
    NumericVector<double>* x_vec = x_system.solution.get();
    NumericVector<double>& X0_vec = x_system.get_vector("INITIAL_COORDINATES");

    NumericVector<double>* x_ghost_vec = x_system.current_local_solution.get();
    x_vec->localize(*x_ghost_vec);
    NumericVector<double>* U_vec = U_system.solution.get();
    NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
    U_vec->localize(*U_ghost_vec);
    const DofMap& dof_map = x_system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);

    std::unique_ptr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<libMesh::Point>& q_point = fe->get_xyz();
    const vector<vector<double> >& phi = fe->get_phi();
    const vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();
    boost::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi;
    dphi_dxi[0] = &fe->get_dphidxi();
    if (NDIM > 2) dphi_dxi[1] = &fe->get_dphideta();

    std::vector<double> U_qp_vec(NDIM);
    std::vector<const std::vector<double>*> var_data(1);
    var_data[0] = &U_qp_vec;
    std::vector<const std::vector<libMesh::VectorValue<double> >*> grad_var_data;
    void* force_fcn_ctx = NULL;

    TensorValue<double> FF_qp;
    boost::multi_array<double, 2> x_node, X_node, U_node;
    VectorValue<double> F_qp, U_qp, x_qp, X_qp, N, n;

    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
        }
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(U_node, *U_ghost_vec, dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(x_qp, qp, x_node, phi);
            jacobian(FF_qp, qp, x_node, dphi);
            interpolate(U_qp, qp, U_node, phi);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_qp_vec[d] = U_qp(d);
            }
            tether_body_force_function_thin(
                F_qp, n, N, FF_qp, x_qp, q_point[qp], elem, 0, var_data, grad_var_data, loop_time, force_fcn_ctx);

            for (int d = 0; d < NDIM; ++d)
            {
                F_integral[d] += F_qp(d) * JxW[qp];
            }
        }
    }
    SAMRAI_MPI::sumReduction(F_integral, NDIM);

    static const double D = 1.0;

    {
        int qp_tot = 0;
        double WSS_L2_norm = 0.0, WSS_max_norm = 0.0;
        double U_L2_norm = 0.0, U_max_norm = 0.0;
        double disp_L2_norm = 0.0, disp_max_norm = 0.0;
        double P_L2_norm = 0.0, P_max_norm = 0.0;
        System& U_system = equation_systems->get_system<System>(IIMethod::VELOCITY_SYSTEM_NAME);
        System& WSS_system = equation_systems->get_system<System>(IIMethod::WSS_IN_SYSTEM_NAME);
        System& P_o_system = equation_systems->get_system<System>(IIMethod::PRESSURE_OUT_SYSTEM_NAME);
        System& P_j_system = equation_systems->get_system<System>(IIMethod::PRESSURE_JUMP_SYSTEM_NAME);

        NumericVector<double>* U_vec = U_system.solution.get();
        NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
        U_vec->localize(*U_ghost_vec);
        DofMap& U_dof_map = U_system.get_dof_map();
        std::vector<std::vector<unsigned int> > U_dof_indices(NDIM);

        NumericVector<double>* WSS_vec = WSS_system.solution.get();
        NumericVector<double>* WSS_ghost_vec = WSS_system.current_local_solution.get();
        WSS_vec->localize(*WSS_ghost_vec);
        DofMap& WSS_dof_map = WSS_system.get_dof_map();
        std::vector<std::vector<unsigned int> > WSS_dof_indices(NDIM);
        std::unique_ptr<FEBase> fe(FEBase::build(dim, WSS_dof_map.variable_type(0)));
        boost::array<VectorValue<double>, 2> dX_dxi, dx_dxi;

        NumericVector<double>* P_o_vec = P_o_system.solution.get();
        NumericVector<double>* P_o_ghost_vec = P_o_system.current_local_solution.get();
        P_o_vec->localize(*P_o_ghost_vec);
        DofMap& P_o_dof_map = P_o_system.get_dof_map();
        std::vector<unsigned int> P_o_dof_indices;

        NumericVector<double>* P_j_vec = P_j_system.solution.get();
        NumericVector<double>* P_j_ghost_vec = P_j_system.current_local_solution.get();
        P_j_vec->localize(*P_j_ghost_vec);
        DofMap& P_j_dof_map = P_j_system.get_dof_map();
        std::vector<unsigned int> P_j_dof_indices;

        VectorValue<double> U_qp, WSS_qp;
        double P_o_qp, P_j_qp, P_i_qp;
        boost::multi_array<double, 2> U_node, WSS_node;
        boost::multi_array<double, 1> P_o_node, P_j_node;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* elem = *el_it;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dof_map.dof_indices(elem, dof_indices[d], d);
                U_dof_map.dof_indices(elem, U_dof_indices[d], d);
                WSS_dof_map.dof_indices(elem, WSS_dof_indices[d], d);
            }
            P_j_dof_map.dof_indices(elem, P_j_dof_indices);
            P_o_dof_map.dof_indices(elem, P_o_dof_indices);
            const int n_qp = qrule->n_points();
            get_values_for_interpolation(U_node, *U_ghost_vec, U_dof_indices);
            get_values_for_interpolation(WSS_node, *WSS_ghost_vec, WSS_dof_indices);
            get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
            get_values_for_interpolation(X_node, X0_vec, dof_indices);
            get_values_for_interpolation(P_j_node, *P_j_ghost_vec, P_j_dof_indices);
            get_values_for_interpolation(P_o_node, *P_o_ghost_vec, P_o_dof_indices);

            for (int qp = 0; qp < n_qp; ++qp)
            {
                interpolate(x_qp, qp, x_node, phi);
                interpolate(X_qp, qp, X_node, phi);
                interpolate(U_qp, qp, U_node, phi);
                interpolate(WSS_qp, qp, WSS_node, phi);
                interpolate(P_o_qp, qp, P_o_node, phi);
                interpolate(P_j_qp, qp, P_j_node, phi);
                P_i_qp = -(P_j_qp - P_o_qp);

                for (unsigned int k = 0; k < NDIM - 1; ++k)
                {
                    interpolate(dx_dxi[k], qp, x_node, *dphi_dxi[k]);
                }
                if (NDIM == 2)
                {
                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                }
                n = (dx_dxi[0].cross(dx_dxi[1])).unit();

                double p_ex_qp = -2. * p_e * x_qp(0) / (L) + p_e;

                double ex_wss[NDIM];
                double ex_U[NDIM];
                ex_wss[0] = -0.5 * p_e * D / L;
                ex_wss[1] = 0.0;
                ex_wss[2] = 0.0;
                ex_U[0] = 0.0;
                ex_U[1] = 0.0;
                ex_U[2] = 0.0;
                if (x_qp(0) > 0.2 * L && x_qp(0) < 0.8 * L)
                {
                    qp_tot += 1;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        U_L2_norm += (U_qp(d) - ex_U[d]) * (U_qp(d) - ex_U[d]) * JxW[qp];
                        U_max_norm = std::max(U_max_norm, std::abs(U_qp(d) - ex_U[d]));
                        WSS_L2_norm += (WSS_qp(d) - ex_wss[d]) * (WSS_qp(d) - ex_wss[d]) * JxW[qp];
                        WSS_max_norm = std::max(WSS_max_norm, std::abs(WSS_qp(d) - ex_wss[d]));
                        disp_L2_norm += (X_qp(d) - x_qp(d)) * (X_qp(d) - x_qp(d)) * JxW[qp];
                        disp_max_norm = std::max(disp_max_norm, std::abs(X_qp(d) - x_qp(d)));
                    }
                    P_L2_norm += std::abs(P_i_qp - p_ex_qp) * std::abs(P_i_qp - p_ex_qp) * JxW[qp];
                    P_max_norm = std::max(P_max_norm, std::abs(P_i_qp - p_ex_qp));
                }
            }
        }
        SAMRAI_MPI::sumReduction(&qp_tot, 1);

        SAMRAI_MPI::sumReduction(&WSS_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&WSS_max_norm, 1);
        SAMRAI_MPI::sumReduction(&U_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&U_max_norm, 1);
        SAMRAI_MPI::sumReduction(&disp_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&disp_max_norm, 1);
        SAMRAI_MPI::sumReduction(&P_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&P_max_norm, 1);

        U_L2_norm = sqrt(U_L2_norm);
        WSS_L2_norm = sqrt(WSS_L2_norm);
        disp_L2_norm = sqrt(disp_L2_norm);
        P_L2_norm = sqrt(P_L2_norm);

        pout << " WSS_L2_norm = " << WSS_L2_norm << "\n\n";
        pout << " WSS_max_norm = " << WSS_max_norm << "\n\n";

        pout << " U_L2_norm = " << U_L2_norm << "\n\n";
        pout << " U_max_norm = " << U_max_norm << "\n\n";

        pout << " disp_L2_norm = " << disp_L2_norm << "\n\n";
        pout << " disp_max_norm = " << disp_max_norm << "\n\n";

        pout << " P_L2_norm = " << P_L2_norm << "\n\n";
        pout << " P_max_norm = " << P_max_norm << "\n\n";
    }

    return;
} // postprocess_data

void
compute_pressure_profile(tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                         const int p_idx,
                         const double /*data_time*/,
                         const string& /*data_dump_dirname*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
    hier_math_ops.resetLevels(finest_ln, finest_ln);
    const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

    const double X_min[3] = { x_loc_min, y_loc_min, z_loc_min };
    const double X_max[3] = { x_loc_max, y_loc_max, z_loc_max };
    // vector<double> pos_values;
    double p_Eulerian_L2_norm = 0.0;
    double p_Eulerian_max_norm = 0.0;

    int qp_tot = 0;
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        tbox::Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            tbox::Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();

            const tbox::Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            const double* const patch_dx = patch_geom->getDx();

            // Entire box containing the required data.
            Box<NDIM> box(IndexUtilities::getCellIndex(
                              &X_min[0], patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper),
                          IndexUtilities::getCellIndex(
                              &X_max[0], patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper));
            // Part of the box on this patch
            Box<NDIM> trim_box = patch_box * box;
            BoxList<NDIM> iterate_box_list = trim_box;

            // Trim the box covered by the finer region
            BoxList<NDIM> covered_boxes;
            if (ln < finest_ln)
            {
                BoxArray<NDIM> refined_region_boxes;
                tbox::Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
                refined_region_boxes = next_finer_level->getBoxes();
                refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM> refined_box = refined_region_boxes[i];
                    const Box<NDIM> covered_box = trim_box * refined_box;
                    covered_boxes.unionBoxes(covered_box);
                }
            }
            iterate_box_list.removeIntersections(covered_boxes);

            // Loop over the boxes and store the location and interpolated value.
            const tbox::Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(p_idx);
            const tbox::Pointer<CellData<NDIM, double> > wgt_cc_data = patch->getPatchData(wgt_cc_idx);

            for (BoxList<NDIM>::Iterator lit(iterate_box_list); lit; lit++)
            {
                const Box<NDIM>& iterate_box = *lit;
                for (Box<NDIM>::Iterator bit(iterate_box); bit; bit++)
                {
                    const CellIndex<NDIM>& cell_idx = *bit;
                    const double p1 = (*p_data)(cell_idx);
                    const double y = patch_x_lower[1] + patch_dx[1] * (cell_idx(1) - patch_lower(1) + 0.5);
                    const double z = patch_x_lower[2] + patch_dx[2] * (cell_idx(2) - patch_lower(2) + 0.5);
                    const double x = patch_x_lower[0] + patch_dx[0] * (cell_idx(0) - patch_lower(0) + 0.5);

                    double p_ex_qp;
                    if (sqrt((y - 0.5 * H) * (y - 0.5 * H) + (z - 0.5 * H) * (z - 0.5 * H)) > D / 2 + patch_dx[1])
                    {
                        p_ex_qp = 0.0;
                    }
                    else if (sqrt((y - 0.5 * H) * (y - 0.5 * H) + (z - 0.5 * H) * (z - 0.5 * H)) < D / 2 - patch_dx[1])
                    {
                        p_ex_qp = -2.0 * p_e * x / L + p_e;
                    }
                    else
                    {
                        p_ex_qp = p1;
                    }

                    if (x > 0.2 * L && x < 0.8 * L)
                    {
                        qp_tot += 1;
                        p_Eulerian_L2_norm +=
                            std::abs(p1 - p_ex_qp) * std::abs(p1 - p_ex_qp) * (*wgt_cc_data)(cell_idx);
                        p_Eulerian_max_norm = std::max(p_Eulerian_max_norm, std::abs(p1 - p_ex_qp));
                    }
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(&qp_tot, 1);
    SAMRAI_MPI::sumReduction(&p_Eulerian_L2_norm, 1);
    SAMRAI_MPI::maxReduction(&p_Eulerian_max_norm, 1);

    p_Eulerian_L2_norm = sqrt(p_Eulerian_L2_norm);

    pout << " p_Eulerian_L2_norm = " << p_Eulerian_L2_norm << "\n\n";
    pout << " p_Eulerian_max_norm = " << p_Eulerian_max_norm << "\n\n";

    return;
} // compute_pressure_profile

void
compute_flow_rate(const tbox::Pointer<PatchHierarchy<NDIM> > hierarchy,
                  const int U_idx,
                  const double loop_time,
                  const int wgt_sc_idx)
{
    vector<double> qsrc;
    qsrc.resize(2);
    std::fill(qsrc.begin(), qsrc.end(), 0.0);
    const double posni[3] = { 0.0, yi, zi };
    const double posno[3] = { L, yo, zo };
    const double rsrc[2] = { 0.5 * D, 0.5 * D };
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        tbox::Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            tbox::Pointer<Patch<NDIM> > patch = level->getPatch(p());
            tbox::Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            if (pgeom->getTouchesRegularBoundary())
            {
                tbox::Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(U_idx);
                tbox::Pointer<SideData<NDIM, double> > wgt_sc_data = patch->getPatchData(wgt_sc_idx);
                const Box<NDIM>& patch_box = patch->getBox();
                const double* const x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                double dV = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    dV *= dx[d];
                }
                double X[NDIM];
                static const int axis = 0;
                for (int side = 0; side <= 1; ++side)
                {
                    const bool is_lower = side == 0;
                    if (pgeom->getTouchesRegularBoundary(axis, side))
                    {
                        // const double rsrc = d_rsrc[side];
                        //
                        VectorNd n;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            n[d] = axis == d ? (is_lower ? -1.0 : +1.0) : 0.0;
                        }
                        Box<NDIM> side_box = patch_box;
                        if (is_lower)
                        {
                            side_box.lower(axis) = patch_box.lower(axis);
                            side_box.upper(axis) = patch_box.lower(axis);
                        }
                        else
                        {
                            side_box.lower(axis) = patch_box.upper(axis) + 1;
                            side_box.upper(axis) = patch_box.upper(axis) + 1;
                        }
                        for (Box<NDIM>::Iterator b(side_box); b; b++)
                        {
                            const hier::Index<NDIM>& i = b();
                            double r_sq = 0.0;
                            for (int d = 0; d < NDIM; ++d)
                            {
                                X[d] =
                                    x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == axis ? 0.0 : 0.5));
                                if (d != axis && side == 0) r_sq += pow(X[d] - posni[d], 2.0);
                                if (d != axis && side == 1) r_sq += pow(X[d] - posno[d], 2.0);
                            }
                            const double r = sqrt(r_sq);
                            if (r <= rsrc[side])
                            {
                                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                {
                                    double dA = n[axis] * dV / dx[axis];
                                    qsrc[side] += (*U_data)(i_s)*dA;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(&qsrc[0], 2);

    if (SAMRAI_MPI::getRank() == 0)
    {
        flow_rate_stream << loop_time << " " << qsrc[0] << " " << qsrc[1] << endl;
    }

    pout << " Qi = " << qsrc[0] << " Qo = " << qsrc[1] << "\n\n";
}

void
postprocess_displacement_data(MeshBase &mesh, System &dX_system)
{
    double max_displacement =0.0;

    NumericVector<double> &dX_vec = *dX_system.solution.get();
    NumericVector<double> &dX_ghost_vec = *dX_system.current_local_solution.get();
    copy_and_synch(dX_vec, dX_ghost_vec);
    DofMap &dX_dof_map = dX_system.get_dof_map();
    std::vector<std::vector<dof_id_type> > dX_dof_indices(NDIM);
    boost::multi_array<double, 2> dX_node;

    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;

        
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dX_dof_map.dof_indices(elem, dX_dof_indices[d], d);
        }


        const int n_basis = static_cast<int>(dX_dof_indices[0].size());
        get_values_for_interpolation(dX_node, dX_ghost_vec, dX_dof_indices);
        for (int k = 0; k < n_basis; ++k)
        {   
	    double current_distance =0.0;
            for (int d = 0; d < NDIM; ++d)
            {
 		current_distance +=std::abs(dX_node[k][d]);  
            }
	    max_displacement = std::max(max_displacement,current_distance);

        }
    }

    SAMRAI_MPI::maxReduction(&max_displacement, 1);
    plog <<"" << max_displacement << std::endl;
    if (SAMRAI_MPI::getRank()==0)
	    max_displacement_stream<< "" <<max_displacement<<std::endl;
} 
