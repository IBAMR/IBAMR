// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Hyperelastic volumetric stabilization utilities following Vadala-Roth,
// Acharya, Patankar, Rossi, and Griffith (2020), Comput. Methods Appl.
// Mech. Engrg. 365:112978, DOI 10.1016/j.cma.2020.112978.

#ifndef included_IBAMR_IBFEVolumetricStabilization
#define included_IBAMR_IBFEVolumetricStabilization

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include <ibtk/libmesh_utilities.h>

#include <tbox/Utilities.h>

#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>

#include <cmath>
#include <string>

namespace IBAMR
{
namespace FEVolumetricStabilization
{
enum class StabilizationFormulation
{
    UNMODIFIED = 0,
    MODIFIED_INVARIANTS = 1,
    DEVIATORIC_PROJECTION = 2
};

inline constexpr int
defaultSpaceDim()
{
#ifdef NDIM
    return NDIM;
#else
    return 3;
#endif
}

inline StabilizationFormulation
stringToFormulation(const std::string& name)
{
    if (name == "UNMODIFIED")
    {
        return StabilizationFormulation::UNMODIFIED;
    }
    if (name == "MODIFIED_INVARIANTS")
    {
        return StabilizationFormulation::MODIFIED_INVARIANTS;
    }
    if (name == "DEVIATORIC_PROJECTION")
    {
        return StabilizationFormulation::DEVIATORIC_PROJECTION;
    }
    TBOX_ERROR("IBAMR::FEVolumetricStabilization::stringToFormulation(): "
               << "unsupported formulation '" << name << "'.\n");
    return StabilizationFormulation::MODIFIED_INVARIANTS;
}

// Evaluate numerical bulk modulus kappa_stab from Eq. (26).
inline double
computeNumericalBulkModulus(const double shear_modulus, const double nu_stab)
{
    if (shear_modulus < 0.0 || nu_stab < -1.0 || nu_stab >= 0.5)
    {
        TBOX_ERROR("IBAMR::FEVolumetricStabilization::computeNumericalBulkModulus(): "
                   << "require shear_modulus >= 0 and -1 <= nu_stab < 0.5.\n");
    }
    return 2.0 * shear_modulus * (1.0 + nu_stab) / (3.0 * (1.0 - 2.0 * nu_stab));
}

inline double
computeMooneyRivlinShearModulus(const double c1, const double c2)
{
    return 2.0 * (c1 + c2);
}

inline double
computePositiveJacobianDeterminant(const libMesh::TensorValue<double>& FF)
{
    const double J = FF.det();
    if (!(J > 0.0))
    {
        TBOX_ERROR("IBAMR::FEVolumetricStabilization::computePositiveJacobianDeterminant(): "
                   << "non-positive Jacobian determinant J = " << J << ".\n");
    }
    return J;
}

// Logarithmic volumetric penalty U(J) = 0.5 * kappa_stab * (ln J)^2, Eq. (25).
inline double
computeLogarithmicPenaltyEnergy(const double J, const double kappa_stab)
{
    if (!(J > 0.0))
    {
        TBOX_ERROR("IBAMR::FEVolumetricStabilization::computeLogarithmicPenaltyEnergy(): "
                   << "require J > 0, got J = " << J << ".\n");
    }
    const double ln_J = std::log(J);
    return 0.5 * kappa_stab * ln_J * ln_J;
}

// PK1 volumetric term kappa_stab * ln(J) * F^{-T}, Eqs. (42), (44), and (45).
inline libMesh::TensorValue<double>
computeLogarithmicPenaltyPK1(const libMesh::TensorValue<double>& FF,
                             const double kappa_stab,
                             const int dim = defaultSpaceDim())
{
    const double J = computePositiveJacobianDeterminant(FF);
    return kappa_stab * std::log(J) * IBTK::tensor_inverse_transpose(FF, dim);
}

// Static pressure P for FEMechanicsBase where PP = -J * P * F^{-T}.
inline double
computeLogarithmicStaticPressure(const double J, const double kappa_stab)
{
    if (!(J > 0.0))
    {
        TBOX_ERROR("IBAMR::FEVolumetricStabilization::computeLogarithmicStaticPressure(): "
                   << "require J > 0, got J = " << J << ".\n");
    }
    return -kappa_stab * std::log(J) / J;
}

// Liu-Hofstetter-Mang energy kappa * (J * ln(J) - J + 1).
inline double
computeLiuHofstetterMangEnergy(const double J, const double kappa_stab)
{
    if (!(J > 0.0))
    {
        TBOX_ERROR("IBAMR::FEVolumetricStabilization::computeLiuHofstetterMangEnergy(): "
                   << "require J > 0, got J = " << J << ".\n");
    }
    return kappa_stab * (J * std::log(J) - J + 1.0);
}

// Liu-Hofstetter-Mang PK1 term kappa * J * ln(J) * F^{-T}.
inline libMesh::TensorValue<double>
computeLiuHofstetterMangPK1(const libMesh::TensorValue<double>& FF,
                            const double kappa_stab,
                            const int dim = defaultSpaceDim())
{
    const double J = computePositiveJacobianDeterminant(FF);
    return kappa_stab * J * std::log(J) * IBTK::tensor_inverse_transpose(FF, dim);
}

// Liu-Hofstetter-Mang static pressure P = -kappa * ln(J).
inline double
computeLiuHofstetterMangStaticPressure(const double J, const double kappa_stab)
{
    if (!(J > 0.0))
    {
        TBOX_ERROR("IBAMR::FEVolumetricStabilization::computeLiuHofstetterMangStaticPressure(): "
                   << "require J > 0, got J = " << J << ".\n");
    }
    return -kappa_stab * std::log(J);
}

// Neo-Hookean elastic energy from Eqs. (41) and (43).
inline double
computeNeoHookeanEnergy(const libMesh::TensorValue<double>& FF,
                        const double shear_modulus,
                        const StabilizationFormulation formulation)
{
    const double J = computePositiveJacobianDeterminant(FF);
    const libMesh::TensorValue<double> CC = FF.transpose() * FF;
    const double I1 = CC.tr();
    switch (formulation)
    {
    case StabilizationFormulation::UNMODIFIED:
        return 0.5 * shear_modulus * (I1 - 3.0);
    case StabilizationFormulation::MODIFIED_INVARIANTS:
    {
        const double I1_bar = std::pow(J, -2.0 / 3.0) * I1;
        return 0.5 * shear_modulus * (I1_bar - 3.0);
    }
    case StabilizationFormulation::DEVIATORIC_PROJECTION:
        TBOX_ERROR("IBAMR::FEVolumetricStabilization::computeNeoHookeanEnergy(): "
                   << "DEVIATORIC_PROJECTION is not hyperelastic.\n");
        return 0.0;
    }
    return 0.0;
}

// Neo-Hookean PK1 stress from Eqs. (42), (44), and (45).
inline libMesh::TensorValue<double>
computeNeoHookeanPK1(const libMesh::TensorValue<double>& FF,
                     const double shear_modulus,
                     const StabilizationFormulation formulation,
                     const int dim = defaultSpaceDim())
{
    const double J = computePositiveJacobianDeterminant(FF);
    if (formulation == StabilizationFormulation::UNMODIFIED)
    {
        return shear_modulus * FF;
    }
    const libMesh::TensorValue<double> CC = FF.transpose() * FF;
    const double I1 = CC.tr();
    const libMesh::TensorValue<double> FF_inv_trans = IBTK::tensor_inverse_transpose(FF, dim);
    const libMesh::TensorValue<double> PP_dev = shear_modulus * (FF - (I1 / 3.0) * FF_inv_trans);
    if (formulation == StabilizationFormulation::MODIFIED_INVARIANTS)
    {
        return std::pow(J, -2.0 / 3.0) * PP_dev;
    }
    return PP_dev;
}

// Mooney-Rivlin elastic energy from Eqs. (46) and (48).
inline double
computeMooneyRivlinEnergy(const libMesh::TensorValue<double>& FF,
                          const double c1,
                          const double c2,
                          const StabilizationFormulation formulation)
{
    const double J = computePositiveJacobianDeterminant(FF);
    const libMesh::TensorValue<double> CC = FF.transpose() * FF;
    const double I1 = CC.tr();
    const double I2 = 0.5 * (I1 * I1 - (CC * CC).tr());
    switch (formulation)
    {
    case StabilizationFormulation::UNMODIFIED:
        return c1 * (I1 - 3.0) + c2 * (I2 - 3.0);
    case StabilizationFormulation::MODIFIED_INVARIANTS:
    {
        const double J_m23 = std::pow(J, -2.0 / 3.0);
        const double I1_bar = J_m23 * I1;
        const double I2_bar = J_m23 * J_m23 * I2;
        return c1 * (I1_bar - 3.0) + c2 * (I2_bar - 3.0);
    }
    case StabilizationFormulation::DEVIATORIC_PROJECTION:
        TBOX_ERROR("IBAMR::FEVolumetricStabilization::computeMooneyRivlinEnergy(): "
                   << "DEVIATORIC_PROJECTION is not hyperelastic.\n");
        return 0.0;
    }
    return 0.0;
}

// Mooney-Rivlin PK1 stress from Eqs. (47), (49), and (50).
inline libMesh::TensorValue<double>
computeMooneyRivlinPK1(const libMesh::TensorValue<double>& FF,
                       const double c1,
                       const double c2,
                       const StabilizationFormulation formulation,
                       const int dim = defaultSpaceDim())
{
    const double J = computePositiveJacobianDeterminant(FF);
    const libMesh::TensorValue<double> CC = FF.transpose() * FF;
    const double I1 = CC.tr();
    const libMesh::TensorValue<double> P1_unmod = 2.0 * c1 * FF;
    const libMesh::TensorValue<double> P2_unmod = 2.0 * c2 * (I1 * FF - FF * CC);
    if (formulation == StabilizationFormulation::UNMODIFIED)
    {
        return P1_unmod + P2_unmod;
    }
    const double I2 = 0.5 * (I1 * I1 - (CC * CC).tr());
    const libMesh::TensorValue<double> FF_inv_trans = IBTK::tensor_inverse_transpose(FF, dim);
    const libMesh::TensorValue<double> P1_dev = P1_unmod - (2.0 * c1 * I1 / 3.0) * FF_inv_trans;
    const libMesh::TensorValue<double> P2_dev = P2_unmod - (4.0 * c2 * I2 / 3.0) * FF_inv_trans;
    if (formulation == StabilizationFormulation::MODIFIED_INVARIANTS)
    {
        const double J_m23 = std::pow(J, -2.0 / 3.0);
        return J_m23 * P1_dev + (J_m23 * J_m23) * P2_dev;
    }
    return P1_dev + P2_dev;
}

inline libMesh::VectorValue<double>
normalizeFiberDirection(const libMesh::VectorValue<double>& fiber_dir)
{
    const double norm = fiber_dir.norm();
    if (!(norm > 0.0))
    {
        TBOX_ERROR("IBAMR::FEVolumetricStabilization::normalizeFiberDirection(): "
                   << "fiber direction vector must be nonzero.\n");
    }
    return fiber_dir / norm;
}

// Modified standard reinforcing energy from Eqs. (53) and (55); I4, I5 stay unmodified.
inline double
computeModifiedReinforcingEnergy(const libMesh::TensorValue<double>& FF,
                                 const double G_T,
                                 const double G_L,
                                 const double E_L,
                                 const libMesh::VectorValue<double>& fiber_dir,
                                 const StabilizationFormulation formulation)
{
    const double W_iso = computeNeoHookeanEnergy(FF, G_T, formulation);
    const libMesh::VectorValue<double> A = normalizeFiberDirection(fiber_dir);
    const libMesh::TensorValue<double> CC = FF.transpose() * FF;
    const libMesh::VectorValue<double> CA = CC * A;
    const double I4 = A * CA;
    const double I5 = CA * CA;
    const double W_aniso =
        0.5 * (G_T - G_L) * (2.0 * I4 - I5 - 1.0) + 0.125 * (E_L + G_T - 4.0 * G_L) * (I4 - 1.0) * (I4 - 1.0);
    return W_iso + W_aniso;
}

// Modified standard reinforcing PK1 stress from Eqs. (54), (56), and (57).
inline libMesh::TensorValue<double>
computeModifiedReinforcingPK1(const libMesh::TensorValue<double>& FF,
                              const double G_T,
                              const double G_L,
                              const double E_L,
                              const libMesh::VectorValue<double>& fiber_dir,
                              const StabilizationFormulation formulation,
                              const int dim = defaultSpaceDim())
{
    const libMesh::TensorValue<double> P_iso = computeNeoHookeanPK1(FF, G_T, formulation, dim);
    const libMesh::VectorValue<double> A = normalizeFiberDirection(fiber_dir);
    libMesh::TensorValue<double> M;
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            M(i, j) = A(i) * A(j);
        }
    }
    const libMesh::TensorValue<double> CC = FF.transpose() * FF;
    const libMesh::VectorValue<double> CA = CC * A;
    const double I4 = A * CA;
    const libMesh::TensorValue<double> FM = FF * M;
    const libMesh::TensorValue<double> P_aniso =
        (G_T - G_L) * (2.0 * FM - FM * CC - FF * CC * M) + 0.5 * (E_L + G_T - 4.0 * G_L) * (I4 - 1.0) * FM;
    return P_iso + P_aniso;
}
} // namespace FEVolumetricStabilization
} // namespace IBAMR

#endif // #ifdef IBAMR_HAVE_LIBMESH
#endif // #ifndef included_IBAMR_IBFEVolumetricStabilization
