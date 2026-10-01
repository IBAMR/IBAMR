// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CartesianCentering.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/NormOps.h>
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/SAMRAIScopedVectorCopy.h>
#include <ibtk/SAMRAIScopedVectorDuplicate.h>
#include <ibtk/muParserCartGridFunction.h>

#include <petscvec.h>

#include <ArrayData.h>
#include <BergerRigoutsos.h>
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <EdgeVariable.h>
#include <FaceVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <NodeVariable.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
// Compare PETSc vector operations applied to wrapped SAMRAI vectors to the
// corresponding SAMRAI vector operations. Coefficients that differ from +/-1 by
// less than sqrt(machine epsilon) must be applied exactly, not as +/-1.
void
check_petsc_vector_ops(SAMRAIVectorReal<NDIM, double>& u_vec, SAMRAIVectorReal<NDIM, double>& f_vec)
{
    SAMRAIScopedVectorDuplicate<double> x_vec(u_vec);
    SAMRAIScopedVectorDuplicate<double> z_vec(f_vec);
    SAMRAIScopedVectorDuplicate<double> y_vec(f_vec);
    SAMRAIScopedVectorDuplicate<double> w_vec(f_vec);
    SAMRAIScopedVectorDuplicate<double> r_vec(f_vec);
    SAMRAIScopedVectorDuplicate<double> e_vec(f_vec);
    Pointer<SAMRAIVectorReal<NDIM, double>> x = x_vec;
    Pointer<SAMRAIVectorReal<NDIM, double>> z = z_vec;
    Pointer<SAMRAIVectorReal<NDIM, double>> y = y_vec;
    Pointer<SAMRAIVectorReal<NDIM, double>> w = w_vec;
    Pointer<SAMRAIVectorReal<NDIM, double>> r = r_vec;
    Pointer<SAMRAIVectorReal<NDIM, double>> e = e_vec;
    x->copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(&u_vec, false), false);
    z->copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false), false);
    const double x_norm = x->maxNorm();
    const double z_norm = z->maxNorm();

    Vec x_petsc = PETScSAMRAIVectorReal::createPETScVector(x);
    Vec z_petsc = PETScSAMRAIVectorReal::createPETScVector(z);
    Vec y_petsc = PETScSAMRAIVectorReal::createPETScVector(y);
    Vec w_petsc = PETScSAMRAIVectorReal::createPETScVector(w);
    std::array<Vec, 2> maxpy_vecs = { x_petsc, z_petsc };

    // Report the error in the PETSc result relative to the operand that is
    // multiplied by the coefficient under test.
    const auto report = [&](const std::string& op,
                            const std::string& coef_name,
                            const Pointer<SAMRAIVectorReal<NDIM, double>>& result,
                            const double operand_norm)
    {
        e->subtract(result, r);
        plog << op << " with coefficient " << coef_name << ": relative error = " << e->maxNorm() / operand_norm << "\n";
    };

    const double delta = 1.0e-8;
    const std::array<std::pair<std::string, double>, 4> coefs = {
        { { "1", 1.0 }, { "-1", -1.0 }, { "1 + delta", 1.0 + delta }, { "-(1 - delta)", -(1.0 - delta) } }
    };
    PetscErrorCode ierr;
    for (const auto& coef : coefs)
    {
        const std::string& name = coef.first;
        const double c = coef.second;
        const std::array<double, 2> maxpy_coefs = { c, c };

        // y = y + c*x
        y->copyVector(z, false);
        r->linearSum(c, x, 1.0, z);
        ierr = VecAXPY(y_petsc, c, x_petsc);
        IBTK_CHKERRQ(ierr);
        report("VecAXPY", name, y, x_norm);

        // y = c*x + 2*y
        y->copyVector(z, false);
        r->linearSum(c, x, 2.0, z);
        ierr = VecAXPBY(y_petsc, c, 2.0, x_petsc);
        IBTK_CHKERRQ(ierr);
        report("VecAXPBY (alpha)", name, y, x_norm);

        // y = 2*x + c*y
        y->copyVector(z, false);
        r->linearSum(2.0, x, c, z);
        ierr = VecAXPBY(y_petsc, 2.0, c, x_petsc);
        IBTK_CHKERRQ(ierr);
        report("VecAXPBY (beta)", name, y, z_norm);

        // y = y + c*x + c*z
        y->copyVector(z, false);
        r->linearSum(c, x, 1.0, z);
        r->axpy(c, z, r);
        ierr = VecMAXPY(y_petsc, 2, maxpy_coefs.data(), maxpy_vecs.data());
        IBTK_CHKERRQ(ierr);
        report("VecMAXPY", name, y, std::max(x_norm, z_norm));

        // y = x + c*y
        y->copyVector(z, false);
        r->linearSum(1.0, x, c, z);
        ierr = VecAYPX(y_petsc, c, x_petsc);
        IBTK_CHKERRQ(ierr);
        report("VecAYPX", name, y, z_norm);

        // w = c*x + y
        y->copyVector(z, false);
        r->linearSum(c, x, 1.0, z);
        ierr = VecWAXPY(w_petsc, c, x_petsc, y_petsc);
        IBTK_CHKERRQ(ierr);
        report("VecWAXPY", name, w, x_norm);
    }

    PETScSAMRAIVectorReal::destroyPETScVector(x_petsc);
    PETScSAMRAIVectorReal::destroyPETScVector(z_petsc);
    PETScSAMRAIVectorReal::destroyPETScVector(y_petsc);
    PETScSAMRAIVectorReal::destroyPETScVector(w_petsc);
}

// A variable of one data centering and depth for vector tests, with the member function of HierarchyMathOps that
// provides its control volume, or nullptr if HierarchyMathOps provides none.
struct TestVariable
{
    std::string label;
    Pointer<SAMRAI::hier::Variable<NDIM>> var;
    int idx;
    int depth;
    int (HierarchyMathOps::*get_weight_index)();
};

// Register a variable with the variable database and return it with its patch data index.
TestVariable
register_test_variable(Pointer<VariableContext> ctx,
                       const std::string& label,
                       Pointer<SAMRAI::hier::Variable<NDIM>> var,
                       const int depth,
                       const int ghost_width,
                       int (HierarchyMathOps::*get_weight_index)())
{
    const int idx =
        VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(var, ctx, IntVector<NDIM>(ghost_width));
    return { label, var, idx, depth, get_weight_index };
}

// Register variables of all five centerings with depths one and two. They are registered before the hierarchy is built
// so that the patch geometry accounts for their ghost widths. Cell and side variables use the control volumes of
// HierarchyMathOps; node, face, and edge variables have none.
std::vector<TestVariable>
register_test_variables(Pointer<VariableContext> ctx)
{
    std::vector<TestVariable> variables;
    for (const auto depth : { 1, 2 })
    {
        const std::string label = " depth " + std::to_string(depth);
        const std::string name = "test_" + std::to_string(depth) + "_";
        variables.push_back(register_test_variable(ctx,
                                                   "cell" + label,
                                                   new CellVariable<NDIM, double>(name + "cell", depth),
                                                   depth,
                                                   1,
                                                   &HierarchyMathOps::getCellWeightPatchDescriptorIndex));
        variables.push_back(register_test_variable(
            ctx, "node" + label, new NodeVariable<NDIM, double>(name + "node", depth), depth, 1, nullptr));
        variables.push_back(register_test_variable(ctx,
                                                   "side" + label,
                                                   new SideVariable<NDIM, double>(name + "side", depth),
                                                   depth,
                                                   1,
                                                   &HierarchyMathOps::getSideWeightPatchDescriptorIndex));
        variables.push_back(register_test_variable(
            ctx, "face" + label, new FaceVariable<NDIM, double>(name + "face", depth), depth, 1, nullptr));
        variables.push_back(register_test_variable(
            ctx, "edge" + label, new EdgeVariable<NDIM, double>(name + "edge", depth), depth, 1, nullptr));
    }
    return variables;
}

// Allocate the patch data of the variables on every level of the hierarchy.
void
allocate_test_variables(Pointer<PatchHierarchy<NDIM>> hierarchy, const std::vector<TestVariable>& variables)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (const auto& variable : variables)
        {
            level->allocatePatchData(variable.idx, 0.0);
        }
    }
}

// Return a vector on the hierarchy whose single component is the variable and its control volume, if it has one.
std::unique_ptr<SAMRAIVectorReal<NDIM, double>>
make_test_vector(Pointer<PatchHierarchy<NDIM>> hierarchy, HierarchyMathOps& hier_math_ops, const TestVariable& variable)
{
    const int cv_idx = variable.get_weight_index ? (hier_math_ops.*variable.get_weight_index)() : invalid_index;
    auto vec = std::make_unique<SAMRAIVectorReal<NDIM, double>>(
        "test " + variable.label, hierarchy, 0, hierarchy->getFinestLevelNumber());
    vec->addComponent(variable.var, variable.idx, cv_idx);
    return vec;
}

// Set the data of component c to smooth functions of the array indices that depend on seed, the direction, and the
// depth, including in the ghost cells. The functions are offset + sin(phase).
template <DataCentering C>
void
fill_component(SAMRAIVectorReal<NDIM, double>& vec, const int c, const double seed, const double offset)
{
    using Traits = CartesianCentering<C>;
    using Data = typename Traits::template Data<double>;
    Pointer<PatchHierarchy<NDIM>> hierarchy = vec.getPatchHierarchy();
    for (int ln = vec.getCoarsestLevelNumber(); ln <= vec.getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Data> data = level->getPatch(p())->getPatchData(vec.getComponentDescriptorIndex(c));
            for (int axis = 0; axis < Traits::num_axes(); ++axis)
            {
                if (!Traits::template has_axis<double>(*data, axis))
                {
                    continue;
                }
                ArrayData<NDIM, double>& array = Traits::template array_data<double>(*data, axis);
                for (Box<NDIM>::Iterator b(array.getBox()); b; b++)
                {
                    for (int depth = 0; depth < array.getDepth(); ++depth)
                    {
                        double phase = seed + 0.3 * axis + 0.7 * depth;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            phase += (0.11 - 0.03 * d) * b()(d);
                        }
                        array(b(), depth) = offset + std::sin(phase);
                    }
                }
            }
        }
    }
}

void
fill_vector(SAMRAIVectorReal<NDIM, double>& vec, const double seed, const double offset = 1.5)
{
    for (int c = 0; c < vec.getNumberOfComponents(); ++c)
    {
        dispatch_data_centering(get_data_centering<double>(*vec.getComponentVariable(c)->getPatchDataFactory()),
                                [&]<DataCentering C>() { fill_component<C>(vec, c, seed + 0.4 * c, offset); });
    }
}

// Return whether a and b differ by a relative difference below tolerance; a difference that is not a number does not.
bool
agree(const double a, const double b, const double tolerance)
{
    const double scale = std::max(std::abs(a), std::abs(b));
    return scale == 0.0 || std::abs(a - b) / scale < tolerance;
}

// Print the L1, L2, and max norms of NormOps for the vector with sorted or plain summation, each line starting with
// prefix, and compare them with those of SAMRAIVectorReal. Return whether they agree.
bool
print_norm_ops(std::ostream& out,
               const std::string& prefix,
               const bool sorted_summation,
               const SAMRAIVectorReal<NDIM, double>& vec)
{
    constexpr double TOLERANCE = 1.0e-12;
    NormOps::setSortedSummation(sorted_summation);
    const double l1_norm = NormOps::L1Norm(&vec);
    const double l2_norm = NormOps::L2Norm(&vec);
    const double max_norm = NormOps::maxNorm(&vec);
    const double expected_l1_norm = vec.L1Norm();
    const double expected_l2_norm = vec.L2Norm();
    const double expected_max_norm = vec.maxNorm();
    const bool agrees = agree(l1_norm, expected_l1_norm, TOLERANCE) && agree(l2_norm, expected_l2_norm, TOLERANCE) &&
                        agree(max_norm, expected_max_norm, TOLERANCE);

    out << "  " << prefix << "L1 norm = " << l1_norm << "\n";
    out << "  " << prefix << "L2 norm = " << l2_norm << "\n";
    out << "  " << prefix << "max norm = " << max_norm << "\n";
    out << "  " << prefix
        << "agree with SAMRAIVectorReal to a relative difference below 1e-12: " << (agrees ? "yes" : "no") << "\n";
    if (!agrees)
    {
        out << "  " << prefix << "SAMRAIVectorReal: L1 norm = " << expected_l1_norm
            << ", L2 norm = " << expected_l2_norm << ", max norm = " << expected_max_norm << "\n";
    }
    return agrees;
}

// Print the norms of NormOps for the vector with sorted summation and then with plain summation, and compare each with
// the norms of SAMRAIVectorReal. Return whether they agree.
bool
check_norm_ops(const std::string& label, const SAMRAIVectorReal<NDIM, double>& vec)
{
    std::ostringstream out;
    out << std::setprecision(12) << label << ":\n";
    const bool sorted_agrees = print_norm_ops(out, "", true, vec);
    const bool plain_agrees = print_norm_ops(out, "plain summation, ", false, vec);
    plog << out.str();
    return sorted_agrees && plain_agrees;
}

// Run check_norm_ops() for a vector of each variable and for a vector with a component of every centering. Report a
// failure only after all comparisons have been printed.
void
check_all_norm_ops(Pointer<PatchHierarchy<NDIM>> hierarchy,
                   HierarchyMathOps& hier_math_ops,
                   const std::vector<TestVariable>& variables)
{
    allocate_test_variables(hierarchy, variables);

    SAMRAIVectorReal<NDIM, double> mixed_vec("test all centerings", hierarchy, 0, hierarchy->getFinestLevelNumber());
    bool all_agree = true;
    for (const auto& variable : variables)
    {
        const auto vec = make_test_vector(hierarchy, hier_math_ops, variable);
        fill_vector(*vec, 0.9 * variable.idx, 0.25);
        all_agree = check_norm_ops(variable.label, *vec) && all_agree;
        if (variable.depth == 1)
        {
            mixed_vec.addComponent(variable.var, variable.idx, vec->getControlVolumeIndex(0));
        }
    }
    all_agree = check_norm_ops("all centerings, depth 1", mixed_vec) && all_agree;
    plog << std::flush;
    if (!all_agree)
    {
        TBOX_ERROR("NormOps norms differ from the norms of SAMRAIVectorReal\n");
    }
}

// Return the maximum of |a[i] - b[i]| for i = 0, ..., n - 1, or infinity if any difference is not a number.
double
max_abs_difference(const double* const a, const double* const b, const int n)
{
    double max_diff = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double diff = std::abs(a[i] - b[i]);
        max_diff = std::max(max_diff, std::isnan(diff) ? std::numeric_limits<double>::infinity() : diff);
    }
    return max_diff;
}

// Compare VecMDot and VecMTDot on wrapped SAMRAI vectors with the sequences of VecDot and VecTDot calls that they fuse;
// the results must be identical. The first vector is a copy of base_vec and the others are copies of operand_vec, which
// differs from base_vec only when the data prevent fusion. The numbers of vectors are below, at, and above the block
// size of four and the group size of eight. Return whether all comparisons are exact.
bool
check_fused_vector_ops(const std::string& label,
                       SAMRAIVectorReal<NDIM, double>& base_vec,
                       SAMRAIVectorReal<NDIM, double>& operand_vec)
{
    constexpr std::array<int, 5> COUNTS = { 1, 3, 4, 5, 9 };
    constexpr int MAX_COUNT = 9;
    std::vector<std::unique_ptr<SAMRAIScopedVectorDuplicate<double>>> duplicates;
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> vecs;
    std::vector<Vec> petsc_vecs;
    for (int k = 0; k <= MAX_COUNT; ++k)
    {
        // Vector 0 is x and vectors 1, ..., MAX_COUNT are the other operands.
        const bool is_operand = k >= 1 && k <= MAX_COUNT;
        duplicates.push_back(std::make_unique<SAMRAIScopedVectorDuplicate<double>>(is_operand ? operand_vec : base_vec,
                                                                                   "fused_" + std::to_string(k)));
        Pointer<SAMRAIVectorReal<NDIM, double>> vec = *duplicates.back();
        fill_vector(*vec, 0.9 * k);
        vecs.push_back(vec);
        petsc_vecs.push_back(PETScSAMRAIVectorReal::createPETScVector(vec));
    }

    std::ostringstream out;
    out << std::setprecision(15) << label << ":\n";
    bool exact = true;
    const auto report = [&](const std::string& name, const std::vector<double>& max_diffs)
    {
        out << "  " << name << " for n = 1, 3, 4, 5, 9:";
        for (const auto max_diff : max_diffs)
        {
            out << ' ' << max_diff;
            exact = exact && max_diff == 0.0;
        }
        out << '\n';
    };

    // VecDot reduces the dot product of each component over the processes before adding the components, which differs
    // in rounding from reducing the sum over the components, as VecMDot does. For a vector with several components, the
    // expected values are therefore the sums of the local dot products of SAMRAI, reduced once.
    const bool single_component = base_vec.getNumberOfComponents() == 1;
    PetscErrorCode ierr;
    std::array<double, MAX_COUNT> dots, expected_dots;
    for (const auto transpose : { false, true })
    {
        std::vector<double> max_diffs;
        for (const auto n : COUNTS)
        {
            ierr = transpose ? VecMTDot(petsc_vecs[0], n, &petsc_vecs[1], dots.data()) :
                               VecMDot(petsc_vecs[0], n, &petsc_vecs[1], dots.data());
            IBTK_CHKERRQ(ierr);
            for (int i = 0; i < n; ++i)
            {
                if (single_component)
                {
                    ierr = transpose ? VecTDot(petsc_vecs[0], petsc_vecs[1 + i], &expected_dots[i]) :
                                       VecDot(petsc_vecs[0], petsc_vecs[1 + i], &expected_dots[i]);
                    IBTK_CHKERRQ(ierr);
                }
                else
                {
                    expected_dots[i] = vecs[0]->dot(vecs[1 + i], /*local_only*/ true);
                }
            }
            if (!single_component)
            {
                IBTK_MPI::sumReduction(expected_dots.data(), n);
            }
            max_diffs.push_back(max_abs_difference(dots.data(), expected_dots.data(), n));
        }
        if (!transpose)
        {
            out << "  VecMDot values for n = 9:";
            for (const auto dot : dots)
            {
                out << ' ' << dot;
            }
            out << '\n';
        }
        report(transpose ? "max |VecMTDot - separate dots|" : "max |VecMDot - separate dots|", max_diffs);
    }

    plog << out.str();
    for (auto& petsc_vec : petsc_vecs)
    {
        PETScSAMRAIVectorReal::destroyPETScVector(petsc_vec);
    }
    return exact;
}

// Run check_fused_vector_ops() for each variable, for a vector with a component of each centering, and for operands
// whose ghost region is wider than that of x, which prevents fusion. Report a failure only after all comparisons have
// been printed.
void
check_all_fused_vector_ops(Pointer<PatchHierarchy<NDIM>> hierarchy,
                           HierarchyMathOps& hier_math_ops,
                           const std::vector<TestVariable>& variables)
{
    allocate_test_variables(hierarchy, variables);

    // The first variable is cell-centered with depth one, and the last is cell-centered with the wider ghost region.
    const TestVariable& wide_variable = variables.back();
    const auto cell_vec = make_test_vector(hierarchy, hier_math_ops, variables.front());
    const auto wide_vec = make_test_vector(hierarchy, hier_math_ops, wide_variable);
    SAMRAIVectorReal<NDIM, double> mixed_vec("fused all centerings", hierarchy, 0, hierarchy->getFinestLevelNumber());
    bool exact = true;
    for (const auto& variable : variables)
    {
        if (&variable == &wide_variable)
        {
            continue;
        }
        const auto vec = make_test_vector(hierarchy, hier_math_ops, variable);
        exact = check_fused_vector_ops(variable.label, *vec, *vec) && exact;
        if (variable.depth == 1)
        {
            mixed_vec.addComponent(variable.var, variable.idx, vec->getControlVolumeIndex(0));
        }
    }
    exact = check_fused_vector_ops("all centerings, depth 1", mixed_vec, mixed_vec) && exact;
    exact = check_fused_vector_ops("cell depth 1 with operands of wider ghost region", *cell_vec, *wide_vec) && exact;
    plog << std::flush;
    if (!exact)
    {
        TBOX_ERROR("fused vector operations differ from the sequences of single-vector operations\n");
    }
}
} // namespace

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_laplace.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        const bool test_copied_vector = input_db->getBoolWithDefault("test_copied_vector", false);
        const bool test_duplicated_vector = input_db->getBoolWithDefault("test_duplicated_vector", false);
        const bool test_standard_vector = !test_copied_vector && !test_duplicated_vector;
        const bool test_petsc_vector_ops = input_db->getBoolWithDefault("test_petsc_vector_ops", false);
        const bool test_norm_ops = input_db->getBoolWithDefault("test_norm_ops", false);
        const bool test_fused_vector_ops = input_db->getBoolWithDefault("test_fused_vector_ops", false);

        // Create major algorithm and data objects that comprise the
        // application. These objects are configured from the input
        // database. Nearly all SAMRAI applications (at least those in IBAMR)
        // start by setting up the same half-dozen objects.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        // We create a variable for every vector we ultimately declare,
        // instead of creating and then cloning vectors. The rationale for
        // this is given below.
        Pointer<CellVariable<NDIM, double>> u_cc_var = new CellVariable<NDIM, double>("u_cc");
        Pointer<CellVariable<NDIM, double>> f_cc_var = new CellVariable<NDIM, double>("f_cc");
        Pointer<CellVariable<NDIM, double>> e_cc_var = new CellVariable<NDIM, double>("e_cc");
        Pointer<CellVariable<NDIM, double>> f_approx_cc_var = new CellVariable<NDIM, double>("f_approx_cc");

        // Internally, SAMRAI keeps track of variables (and their
        // corresponding vectors, data, etc.) by converting them to
        // indices. Here we get the indices after notifying the variable
        // database about them.
        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(1));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector<NDIM>(1));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(1));
        const int f_approx_cc_idx = var_db->registerVariableAndContext(f_approx_cc_var, ctx, IntVector<NDIM>(1));

        std::vector<TestVariable> test_variables;
        if (test_norm_ops || test_fused_vector_ops)
        {
            test_variables = register_test_variables(ctx);
        }
        if (test_fused_vector_ops)
        {
            test_variables.push_back(register_test_variable(ctx,
                                                            "wide cell",
                                                            new CellVariable<NDIM, double>("test_wide_cell"),
                                                            1,
                                                            2,
                                                            &HierarchyMathOps::getCellWeightPatchDescriptorIndex));
        }

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while ((gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }

        const int finest_level = patch_hierarchy->getFinestLevelNumber();

        // Allocate data for each variable on each level of the patch
        // hierarchy.
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
            level->allocatePatchData(f_approx_cc_idx, 0.0);
        }

        // By default, the norms defined on SAMRAI vectors are vectors in R^n:
        // however, in IBAMR we almost always want to use a norm that
        // corresponds to a numerical quadrature. To do this we have to
        // associate each vector with a set of cell-centered volumes. Rather
        // than set this up manually, we rely on an IBTK utility class that
        // computes this (as well as many other things!). These values are
        // known as `cell weights' in this context, so we get the ID of the
        // associated data by asking for that. Behind the scenes
        // HierarchyMathOps sets up the necessary cell-centered variables and
        // registers them with the usual SAMRAI objects: all we need to do is
        // ask for the ID. Due to the way SAMRAI works these calls must occur
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int cv_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        if (test_norm_ops)
        {
            check_all_norm_ops(patch_hierarchy, hier_math_ops, test_variables);
            return EXIT_SUCCESS;
        }
        if (test_fused_vector_ops)
        {
            check_all_fused_vector_ops(patch_hierarchy, hier_math_ops, test_variables);
            return EXIT_SUCCESS;
        }

        // SAMRAI patches do not store data as a single contiguous arrays;
        // instead, each hierarchy contains several contiguous arrays. Hence,
        // to do linear algebra, we rely on SAMRAI's own vector class which
        // understands these relationships. We begin by initializing each
        // vector with the patch hierarchy:
        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, finest_level);
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, finest_level);
        SAMRAIVectorReal<NDIM, double> f_standard("f_approx", patch_hierarchy, 0, finest_level);

        f_vec.addComponent(f_cc_var, f_cc_idx, cv_cc_idx);
        SAMRAIScopedVectorDuplicate<double> f_duplicated(f_vec);
        SAMRAIScopedVectorCopy<double> f_copied(f_vec);
        SAMRAIVectorReal<NDIM, double>* f_approx_vec_ptr = nullptr;

        if (test_copied_vector)
        {
            f_approx_vec_ptr = &static_cast<SAMRAIVectorReal<NDIM, double>&>(f_copied);
        }
        else if (test_duplicated_vector)
        {
            f_approx_vec_ptr = &static_cast<SAMRAIVectorReal<NDIM, double>&>(f_duplicated);
        }
        else if (test_standard_vector)
        {
            f_standard.addComponent(f_approx_cc_var, f_approx_cc_idx, cv_cc_idx);
            f_approx_vec_ptr = &f_standard;
        }
        else
        {
            TBOX_ERROR("unknown test configuration - should be copied, duplicated, or standard");
        }

        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, finest_level);

        u_vec.addComponent(u_cc_var, u_cc_idx, cv_cc_idx);
        e_vec.addComponent(e_cc_var, e_cc_idx, cv_cc_idx);

        u_vec.setToScalar(0.0, false);
        f_vec.setToScalar(0.0, false);
        if (!test_duplicated_vector)
        {
            f_approx_vec_ptr->setToScalar(0.0, false);
        }
        e_vec.setToScalar(0.0, false);

        // Next, we use functions defined with muParser to set up the right
        // hand side and solution. These functions are read from the input
        // database and can be changed without recompiling.
        {
            muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
            muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);

            u_fcn.setDataOnPatchHierarchy(u_cc_idx, u_cc_var, patch_hierarchy, 0.0);
            f_fcn.setDataOnPatchHierarchy(f_cc_idx, f_cc_var, patch_hierarchy, 0.0);
        }

        if (test_petsc_vector_ops) check_petsc_vector_ops(u_vec, f_vec);

        // Compute -L*u = f.
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCConstant(0.0);
        poisson_spec.setDConstant(-1.0);
        RobinBcCoefStrategy<NDIM>* bc_coef = nullptr;
        CCLaplaceOperator laplace_op("laplace op");
        laplace_op.setPoissonSpecifications(poisson_spec);
        laplace_op.setPhysicalBcCoef(bc_coef);
        laplace_op.initializeOperatorState(u_vec, f_vec);
        if (test_copied_vector)
        {
            laplace_op.apply(u_vec, f_copied);
        }
        else if (test_duplicated_vector)
        {
            laplace_op.apply(u_vec, f_duplicated);
        }
        else
        {
            laplace_op.apply(u_vec, f_standard);
        }

        // Compute error and print error norms. Here we create temporary smart
        // pointers that will not delete the underlying object since the
        // second argument to the constructor is false.
        if (test_copied_vector)
        {
            e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false), f_copied);
        }
        else if (test_duplicated_vector)
        {
            e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false), f_duplicated);
        }
        else
        {
            e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false),
                           Pointer<SAMRAIVectorReal<NDIM, double>>(&f_standard, false));
        }
        const double max_norm = e_vec.maxNorm();
        const double l2_norm = e_vec.L2Norm();
        const double l1_norm = e_vec.L1Norm();

        plog << "|e|_oo = " << max_norm << "\n";
        plog << "|e|_2  = " << l2_norm << "\n";
        plog << "|e|_1  = " << l1_norm << "\n";

        {
            std::ostringstream out;
            for (int ln = 0; ln <= finest_level; ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM>> patch_level = patch_hierarchy->getPatchLevel(ln);
                out << std::setprecision(20);
                out << "rank: " << IBTK_MPI::getRank() << " level: " << ln << " boxes:\n";
                for (typename hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
                {
                    const hier::Box<NDIM> box = patch_level->getPatch(p())->getBox();
                    Pointer<CartesianPatchGeometry<NDIM>> patch_geometry =
                        patch_level->getPatch(p())->getPatchGeometry();
                    out << "  " << box << '\n';

                    out << "  x_lo = ";
                    for (int d = 0; d < NDIM - 1; ++d) out << patch_geometry->getXLower()[d] << ", ";
                    out << patch_geometry->getXLower()[NDIM - 1] << std::endl;

                    out << "  x_up = ";
                    for (int d = 0; d < NDIM - 1; ++d) out << patch_geometry->getXUpper()[d] << ", ";
                    out << patch_geometry->getXUpper()[NDIM - 1] << std::endl;

                    out << "  dx   = ";
                    for (int d = 0; d < NDIM - 1; ++d) out << patch_geometry->getDx()[d] << ", ";
                    out << patch_geometry->getDx()[NDIM - 1] << std::endl;
                }

                if (ln != finest_level) out << std::endl;
            }
            print_strings_on_plog_0(out.str());
        }

        // Finally, we clean up the output by setting error values on patches
        // on coarser levels which are covered by finer levels to zero.
        for (int ln = 0; ln < finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM>> next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            BoxArray<NDIM> refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                const Patch<NDIM>& patch = *level->getPatch(p());
                const Box<NDIM>& patch_box = patch.getBox();
                Pointer<CellData<NDIM, double>> e_cc_data = patch.getPatchData(e_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM>& refined_box = refined_region_boxes[i];
                    // Box::operator* returns the intersection of two boxes.
                    const Box<NDIM>& intersection = patch_box * refined_box;
                    if (!intersection.empty())
                    {
                        e_cc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }

        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "SCALAR", u_cc_idx);
        visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "SCALAR", f_cc_idx);
        visit_data_writer->registerPlotQuantity(f_approx_cc_var->getName(), "SCALAR", f_approx_cc_idx);
        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "SCALAR", e_cc_idx);
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
    }
} // run_example
