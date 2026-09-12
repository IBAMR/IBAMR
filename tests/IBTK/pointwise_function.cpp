// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2026 by the IBAMR developers
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
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/CartGridPointwiseFunction.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

#include <tbox/Logger.h>
#include <tbox/MemoryDatabase.h>

#include <BergerRigoutsos.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <ibtk/app_namespaces.h>

namespace
{
constexpr double INITIAL_TIME = 0.5;
constexpr double TRANSFORM_TIME = 1.25;
using Functions = std::array<Pointer<CartGridFunction>, 2>;

double
coordinate_value(const VectorNd& x, const double time)
{
    double value = time;
    for (int d = 0; d < NDIM; ++d)
    {
        value += (d + 1) * x[d];
    }
    return value;
}

int
axis_number(const int axis)
{
    return axis == invalid_index ? 0 : axis + 1;
}

std::string
coordinate_expression(const std::string& time)
{
    std::string expression = time;
    for (int d = 0; d < NDIM; ++d)
    {
        expression += "+" + std::to_string(d + 1) + "*X_" + std::to_string(d);
    }
    return "(" + expression + ")";
}

Functions
make_scalar_functions(Pointer<Variable<NDIM>> var)
{
    return { make_cart_grid_pointwise_function<double>(
                 "scalar initialization",
                 var,
                 [offset = std::make_unique<double>(0.0)](const VectorNd& x, double t, int d, int axis)
                 { return coordinate_value(x, t) + 10 * d + 100 * axis_number(axis) + *offset; }),
             make_cart_grid_pointwise_function<double>("scalar transformation",
                                                       var,
                                                       [](double q, const VectorNd& x, double t, int d, int axis) {
                                                           return 2 * q + coordinate_value(x, t) + d +
                                                                  axis_number(axis);
                                                       }) };
}

template <typename Value>
Functions
make_vector_functions(Pointer<Variable<NDIM>> var, const int depth)
{
    return { make_cart_grid_pointwise_function<Value>(
                 "vector initialization",
                 var,
                 [depth](const VectorNd& x, double t, int d, int axis) -> Value
                 {
                     if (d != 0)
                     {
                         TBOX_ERROR("Whole-vector callbacks must receive depth zero\n");
                     }
                     Value q;
                     if constexpr (std::is_same_v<Value, VectorXd>)
                     {
                         q.resize(depth);
                     }
                     for (int k = 0; k < depth; ++k)
                     {
                         q[k] = coordinate_value(x, t) + 10 * k + 100 * axis_number(axis);
                     }
                     return q;
                 }),
             make_cart_grid_pointwise_function<Value>("vector transformation",
                                                      var,
                                                      [](const Value& q, const VectorNd&, double, int, int)
                                                      { return q.reverse() + 2.0 * q; }) };
}

MatrixNd
base_tensor(const bool symmetric)
{
    MatrixNd q;
    for (int i = 0; i < NDIM; ++i)
    {
        for (int j = 0; j < NDIM; ++j)
        {
            q(i, j) = symmetric ? (i == j ? 5.0 + i : 0.25 * (i + j + 1)) : 3.0 * i + j + 1.0;
        }
    }
    return q;
}

Functions
make_tensor_functions(Pointer<Variable<NDIM>> var, const TensorStorage storage)
{
    return { make_cart_grid_pointwise_function<MatrixNd>(
                 "tensor initialization",
                 var,
                 [base = base_tensor(storage == TensorStorage::SYMMETRIC)](
                     const VectorNd& x, double t, int d, int axis) -> MatrixNd
                 {
                     if (d != 0)
                     {
                         TBOX_ERROR("Whole-tensor callbacks must receive depth zero\n");
                     }
                     return base.array() + coordinate_value(x, t) + 100 * axis_number(axis);
                 },
                 storage),
             make_cart_grid_pointwise_function<MatrixNd>(
                 "tensor transformation",
                 var,
                 [](const MatrixNd& q, const VectorNd&, double, int, int) { return q * q; },
                 storage) };
}

// Independent component formulas for the muParser reference. Tensor packing is
// specified explicitly here instead of calling the implementation's helpers.
std::vector<std::string>
reference_expressions(const std::string& kind, const int depth, const bool staggered, const bool transformed)
{
    std::vector<std::string> expressions;
    const int n_axes = staggered ? NDIM : 1;
    for (int d = 0; d < depth; ++d)
    {
        for (int axis = 0; axis < n_axes; ++axis)
        {
            const int a = staggered ? axis + 1 : 0;
            const std::string coordinate = coordinate_expression(transformed ? "0.5" : "t");
            const auto component = [&](int k)
            { return "(" + coordinate + "+" + std::to_string(10 * k + 100 * a) + ")"; };
            std::string expression;
            if (kind == "scalar")
            {
                expression = transformed ?
                                 "2*" + component(d) + "+" + coordinate_expression("t") + "+" + std::to_string(d + a) :
                                 component(d);
            }
            else if (kind == "vector" || kind == "general")
            {
                expression = transformed ? component(depth - 1 - d) + "+2*" + component(d) : component(d);
            }
            else
            {
                const bool symmetric = kind == "symmetric tensor";
                const MatrixNd base = base_tensor(symmetric);
                int row = d / NDIM;
                int col = d % NDIM;
                if (symmetric)
                {
#if NDIM == 2
                    const int rows[] = { 0, 1, 0 };
                    const int cols[] = { 0, 1, 1 };
#else
                    const int rows[] = { 0, 1, 2, 1, 0, 0 };
                    const int cols[] = { 0, 1, 2, 2, 2, 1 };
#endif
                    row = rows[d];
                    col = cols[d];
                }
                const auto tensor_component = [&](int i, int j)
                { return "(" + coordinate + "+" + std::to_string(base(i, j) + 100 * a) + ")"; };
                if (transformed)
                {
                    expression = "0";
                    for (int k = 0; k < NDIM; ++k)
                    {
                        expression += "+" + tensor_component(row, k) + "*" + tensor_component(k, col);
                    }
                }
                else
                {
                    expression = tensor_component(row, col);
                }
            }
            expressions.push_back(expression);
        }
    }
    return expressions;
}

void
set_reference(Pointer<PatchHierarchy<NDIM>> hierarchy,
              const int data_idx,
              Pointer<Variable<NDIM>> var,
              const std::vector<std::string>& expressions,
              const double time)
{
    Pointer<Database> db = new MemoryDatabase("reference");
    for (unsigned int d = 0; d < expressions.size(); ++d)
    {
        db->putString("function_" + std::to_string(d), expressions[d]);
    }
    muParserCartGridFunction reference("reference", db, hierarchy->getGridGeometry());
    reference.setDataOnPatchHierarchy(data_idx, var, hierarchy, time);
}

double
array_error(const ArrayData<NDIM, double>& data, const ArrayData<NDIM, double>& reference)
{
    double error = 0.0;
    for (Box<NDIM>::Iterator it(data.getBox()); it; it++)
    {
        for (int d = 0; d < data.getDepth(); ++d)
        {
            const double actual = data(it(), d);
            const double expected = reference(it(), d);
            if (std::isnan(expected))
            {
                TBOX_ASSERT(std::isnan(actual));
            }
            else
            {
                TBOX_ASSERT(std::isfinite(actual));
                error = std::max(error, std::abs(actual - expected));
            }
        }
    }
    return error;
}

template <typename Data>
double
compute_error(Pointer<PatchHierarchy<NDIM>> hierarchy, const int data_idx, const int reference_idx)
{
    double error = 0.0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator it(level); it; it++)
        {
            Pointer<Data> data = level->getPatch(it())->getPatchData(data_idx);
            Pointer<Data> reference = level->getPatch(it())->getPatchData(reference_idx);
            if constexpr (std::is_same_v<Data, CellData<NDIM, double>> || std::is_same_v<Data, NodeData<NDIM, double>>)
            {
                error = std::max(error, array_error(data->getArrayData(), reference->getArrayData()));
            }
            else
            {
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    if constexpr (std::is_same_v<Data, SideData<NDIM, double>>)
                    {
                        if (!data->getDirectionVector()(axis))
                        {
                            continue;
                        }
                    }
                    error = std::max(error, array_error(data->getArrayData(axis), reference->getArrayData(axis)));
                }
            }
        }
    }
    return IBTK_MPI::maxReduction(error);
}

template <typename Data, typename Var>
void
run_case(Pointer<PatchHierarchy<NDIM>> hierarchy,
         const std::string& centering,
         const std::string& kind,
         const int depth,
         const bool partial = false)
{
    const std::string name = centering + " " + kind;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("pointwise");
    // Select centering from another variable; the actual data supply depth and directions.
    Pointer<Var> selector_var = new Var(name + " selector");
    Functions functions;
    if (kind == "scalar")
    {
        functions = make_scalar_functions(selector_var);
    }
    else if (kind == "vector")
    {
        functions = make_vector_functions<VectorNd>(selector_var, depth);
    }
    else if (kind == "general")
    {
        functions = make_vector_functions<VectorXd>(selector_var, depth);
    }
    else
    {
        functions =
            make_tensor_functions(selector_var, kind == "full tensor" ? TensorStorage::FULL : TensorStorage::SYMMETRIC);
    }
    selector_var.setNull();
    Pointer<Var> var;
    if constexpr (std::is_same_v<Var, SideVariable<NDIM, double>>)
    {
        var = new Var(name, depth, true, partial ? NDIM - 1 : -1);
    }
    else
    {
        var = new Var(name, depth);
    }
    Pointer<Var> reference_var = new Var(name + " reference", depth);
    IntVector<NDIM> ghosts(1);
    ghosts(1) = 2;
    const int data_idx = var_db->registerVariableAndContext(var, context, ghosts);
    const int reference_idx = var_db->registerVariableAndContext(reference_var, context, ghosts);
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->allocatePatchData(data_idx);
        level->allocatePatchData(reference_idx);
        for (PatchLevel<NDIM>::Iterator it(level); it; it++)
        {
            Pointer<Data> data = level->getPatch(it())->getPatchData(data_idx);
            Pointer<Data> reference = level->getPatch(it())->getPatchData(reference_idx);
            data->fillAll(std::numeric_limits<double>::quiet_NaN());
            reference->fillAll(std::numeric_limits<double>::quiet_NaN());
        }
    }
    const bool staggered =
        !std::is_same_v<Data, CellData<NDIM, double>> && !std::is_same_v<Data, NodeData<NDIM, double>>;
    set_reference(
        hierarchy, reference_idx, reference_var, reference_expressions(kind, depth, staggered, false), INITIAL_TIME);
    functions[0]->setDataOnPatchHierarchy(data_idx, var, hierarchy, INITIAL_TIME);
    const double initialization_error = compute_error<Data>(hierarchy, data_idx, reference_idx);
    functions[1]->setDataOnPatchHierarchy(data_idx, var, hierarchy, TRANSFORM_TIME);
    set_reference(
        hierarchy, reference_idx, reference_var, reference_expressions(kind, depth, staggered, true), TRANSFORM_TIME);
    const double transformation_error = compute_error<Data>(hierarchy, data_idx, reference_idx);

    CartGridFunctionSet sum(name + " sum");
    sum.addFunction(functions[0]);
    sum.addFunction(functions[0]);
    std::vector<std::string> sum_expressions = reference_expressions(kind, depth, staggered, false);
    for (std::string& expression : sum_expressions)
    {
        expression = "2*(" + expression + ")";
    }
    set_reference(hierarchy, reference_idx, reference_var, sum_expressions, INITIAL_TIME);
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator it(level); it; it++)
        {
            sum.setDataOnPatch(data_idx, var, level->getPatch(it()), INITIAL_TIME, true, level);
        }
    }
    const double patch_sum_error = compute_error<Data>(hierarchy, data_idx, reference_idx);
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        sum.setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(ln), INITIAL_TIME, true);
    }
    const double level_sum_error = compute_error<Data>(hierarchy, data_idx, reference_idx);
    plog << name << " errors (initialize, transform, patch sum, level sum): " << initialization_error << ' '
         << transformation_error << ' ' << patch_sum_error << ' ' << level_sum_error << '\n';
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(data_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(reference_idx);
    }
    var_db->removePatchDataIndex(data_idx);
    var_db->removePatchDataIndex(reference_idx);
}

// Preserve the actual abort diagnostic while omitting source paths and line numbers.
class ErrorAppender : public Logger::Appender
{
public:
    /*! \brief Flush the abort diagnostic directly to the test output before termination. */
    void logMessage(const std::string& message, const std::string&, const int) override
    {
        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream output("output");
            output << message.c_str() << std::flush;
        }
    }
};

void
run_error_case(Pointer<PatchHierarchy<NDIM>> hierarchy, const std::string& error)
{
    Logger::getInstance()->setAbortAppender(new ErrorAppender());
    Pointer<Variable<NDIM>> selector_var = new CellVariable<NDIM, double>("selector");
    Pointer<CartGridFunction> function;
    const int depth = NDIM;
    if (error == "dynamic_shape")
    {
        function = make_cart_grid_pointwise_function<VectorXd>(
            "bad shape", selector_var, [](const VectorNd&, double, int, int) { return VectorXd::Zero(NDIM + 1); });
    }
    else if (error == "centering" || error == "data_type")
    {
        function = make_scalar_functions(selector_var)[0];
    }
    else
    {
        TBOX_ERROR("Unknown error case\n");
    }
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM>> var;
    if (error == "data_type")
    {
        var = new CellVariable<NDIM, int>("invalid", depth);
    }
    else if (error == "centering")
    {
        var = new SideVariable<NDIM, double>("invalid", depth);
    }
    else
    {
        var = new CellVariable<NDIM, double>("invalid", depth);
    }
    const int idx = var_db->registerVariableAndContext(var, var_db->getContext("invalid"), IntVector<NDIM>(0));
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    level->allocatePatchData(idx);
    function->setDataOnPatchHierarchy(idx, var, hierarchy, INITIAL_TIME);
    // A zero exit status must fail an expect_error test if no diagnostic occurred.
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv);
    Logger::getInstance()->setWarning(false);
    Pointer<CartesianGridGeometry<NDIM>> geometry =
        new CartesianGridGeometry<NDIM>("geometry", app->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("hierarchy", geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector =
        new StandardTagAndInitialize<NDIM>("tagging", nullptr, app->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load = new LoadBalancer<NDIM>("load", app->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
        "gridding", app->getComponentDatabase("GriddingAlgorithm"), error_detector, boxes, load);
    // Register the maximum test ghost width before constructing the hierarchy.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellVariable<NDIM, double>> width_var = new CellVariable<NDIM, double>("ghost width");
    var_db->registerVariableAndContext(width_var, var_db->getContext("width"), IntVector<NDIM>(2));
    gridding->makeCoarsestLevel(hierarchy, 0.0);
    int ln = 0;
    while (gridding->levelCanBeRefined(ln))
    {
        gridding->makeFinerLevel(hierarchy, 0.0, 0.0, 1);
        if (!hierarchy->finerLevelExists(ln))
        {
            break;
        }
        ++ln;
    }
    const std::string error = app->getInputDatabase()->getStringWithDefault("error_case", "");
    if (!error.empty())
    {
        run_error_case(hierarchy, error);
        return 0;
    }
    if (hierarchy->getFinestLevelNumber() != 1)
    {
        TBOX_ERROR("This test requires two hierarchy levels\n");
    }
    plog << std::scientific << std::setprecision(12);
    // Exercise each centering's geometry with scalars, then vector and tensor
    // packing on representative layouts.
    run_case<CellData<NDIM, double>, CellVariable<NDIM, double>>(hierarchy, "cell", "scalar", NDIM + 2);
    run_case<NodeData<NDIM, double>, NodeVariable<NDIM, double>>(hierarchy, "node", "scalar", NDIM + 2);
    run_case<SideData<NDIM, double>, SideVariable<NDIM, double>>(hierarchy, "side", "scalar", NDIM + 2);
    run_case<FaceData<NDIM, double>, FaceVariable<NDIM, double>>(hierarchy, "face", "scalar", NDIM + 2);
    run_case<EdgeData<NDIM, double>, EdgeVariable<NDIM, double>>(hierarchy, "edge", "scalar", NDIM + 2);
    run_case<SideData<NDIM, double>, SideVariable<NDIM, double>>(hierarchy, "partial side", "scalar", NDIM + 2, true);
    run_case<NodeData<NDIM, double>, NodeVariable<NDIM, double>>(hierarchy, "node", "vector", NDIM);
    run_case<FaceData<NDIM, double>, FaceVariable<NDIM, double>>(hierarchy, "face", "general", 2 * NDIM + 1);
    run_case<CellData<NDIM, double>, CellVariable<NDIM, double>>(hierarchy, "cell", "full tensor", NDIM * NDIM);
    run_case<SideData<NDIM, double>, SideVariable<NDIM, double>>(
        hierarchy, "side", "symmetric tensor", NDIM * (NDIM + 1) / 2);
    return 0;
}
