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
#include <ibtk/CartGridPatchwiseFunction.h>
#include <ibtk/CartGridPointwiseFunction.h>
#include <ibtk/CartesianCentering.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

#include <tbox/Logger.h>
#include <tbox/MemoryDatabase.h>

#include <BergerRigoutsos.h>
#include <CartesianPatchGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <OutersideVariable.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
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

class ScalarResult
{
public:
    /*! \brief Store the scalar callback result. */
    explicit ScalarResult(const double value) : d_value(value)
    {
    }

    /*! \brief Require conversion from the returned value rather than an lvalue. */
    operator double() &&
    {
        return d_value;
    }

private:
    double d_value;
};

Functions
make_scalar_functions(Pointer<Variable<NDIM>> var)
{
    // Copy a const source functor and preserve its rvalue-only scalar conversion.
    const auto transform = [result = 0.0](double q, const VectorNd& x, double t, int d, int axis) mutable
    {
        result = 2 * q + coordinate_value(x, t) + d + axis_number(axis);
        return ScalarResult(result);
    };
    return { make_cart_grid_pointwise_function<double>(
                 "scalar initialization",
                 var,
                 [offset = std::make_unique<double>(0.0)](const VectorNd& x, double t, int d, int axis)
                 { return coordinate_value(x, t) + 10 * d + 100 * axis_number(axis) + *offset; }),
             make_cart_grid_pointwise_function<double>("scalar transformation", var, transform) };
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
    // Exercise const lvalue copying in the tensor factory as well.
    const auto initialize = [base = base_tensor(storage == TensorStorage::SYMMETRIC)](
                                const VectorNd& x, double t, int d, int axis) mutable -> MatrixNd
    {
        if (d != 0)
        {
            TBOX_ERROR("Whole-tensor callbacks must receive depth zero\n");
        }
        return base.array() + coordinate_value(x, t) + 100 * axis_number(axis);
    };
    return { make_cart_grid_pointwise_function<MatrixNd>("tensor initialization", var, initialize, storage),
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

// Reference geometry box for a centering, computed directly with SAMRAI.
template <DataCentering C>
Box<NDIM>
geometry_box(const Box<NDIM>& cell_box, const int axis)
{
    if constexpr (C == DataCentering::CELL)
    {
        return cell_box;
    }
    else if constexpr (C == DataCentering::NODE)
    {
        return NodeGeometry<NDIM>::toNodeBox(cell_box);
    }
    else if constexpr (C == DataCentering::SIDE)
    {
        return SideGeometry<NDIM>::toSideBox(cell_box, axis);
    }
    else if constexpr (C == DataCentering::FACE)
    {
        return FaceGeometry<NDIM>::toFaceBox(cell_box, axis);
    }
    else
    {
        return EdgeGeometry<NDIM>::toEdgeBox(cell_box, axis);
    }
}

// Reference storage array for a centering, taken directly from the SAMRAI accessor.
template <DataCentering C>
const ArrayData<NDIM, double>&
direct_array_data(const typename CartesianCentering<C>::template Data<double>& data, const int axis)
{
    if constexpr (CartesianCentering<C>::is_staggered())
    {
        return data.getArrayData(axis);
    }
    else
    {
        return data.getArrayData();
    }
}

void
print_box(const Box<NDIM>& box)
{
    plog << "lower " << box.lower() << " upper " << box.upper();
}

template <DataCentering C>
void
report_centering(const Box<NDIM>& cell_box, const int depth)
{
    using Layout = CartesianCentering<C>;
    using Data = typename Layout::template Data<double>;
    Data data(cell_box, depth, IntVector<NDIM>(0));
    const Data& const_data = data;
    plog << enum_to_string(C) << ": staggered " << Layout::is_staggered() << ", num_axes " << Layout::num_axes()
         << '\n';
    for (int axis = 0; axis < Layout::num_axes(); ++axis)
    {
        const Box<NDIM> box = Layout::index_box(cell_box, axis);
        const ArrayData<NDIM, double>& array = Layout::template array_data<double>(data, axis);
        const ArrayData<NDIM, double>& const_array = Layout::template array_data<double>(const_data, axis);
        const ArrayData<NDIM, double>& direct = direct_array_data<C>(const_data, axis);
        plog << "  axis " << axis << ": index_box ";
        print_box(box);
        plog << "; equals SAMRAI geometry box " << (box == geometry_box<C>(cell_box, axis)) << "; array_data box ";
        print_box(array.getBox());
        plog << " depth " << array.getDepth() << ", equals index_box " << (array.getBox() == box)
             << "; same array as SAMRAI accessor (non-const, const) " << (&array == &direct) << ' '
             << (&const_array == &direct) << '\n';
    }
}

template <typename T>
void
report_factory(const std::string& label, const int depth, Pointer<Variable<NDIM>> var)
{
    const SAMRAI::hier::PatchDataFactory<NDIM>& factory = *var->getPatchDataFactory();
    const std::optional<DataCentering> centering = find_data_centering<T>(factory);
    plog << label << " (depth " << depth << "): find_data_centering has value " << centering.has_value();
    if (centering)
    {
        plog << ", centering " << enum_to_string(*centering) << " (" << static_cast<int>(*centering)
             << "), get_data_centering " << enum_to_string(get_data_centering<T>(factory));
    }
    plog << '\n';
}

void
run_centering_helpers()
{
    // A small cell box that is not a cube, so that every axis gives a distinct box.
    SAMRAI::hier::Index<NDIM> lower;
    SAMRAI::hier::Index<NDIM> upper;
    for (int d = 0; d < NDIM; ++d)
    {
        lower(d) = d + 1;
        upper(d) = 2 * d + 3;
    }
    const Box<NDIM> cell_box(lower, upper);
    const int depth = 2;
    plog << "cell box ";
    print_box(cell_box);
    plog << '\n';
    for (const DataCentering centering :
         { DataCentering::CELL, DataCentering::NODE, DataCentering::SIDE, DataCentering::FACE, DataCentering::EDGE })
    {
        dispatch_data_centering(centering, [&]<DataCentering C>() { report_centering<C>(cell_box, depth); });
    }
    report_factory<double>("cell double", depth, new CellVariable<NDIM, double>("cell", depth));
    report_factory<double>("node double", depth, new NodeVariable<NDIM, double>("node", depth));
    report_factory<double>("side double", depth, new SideVariable<NDIM, double>("side", depth));
    report_factory<double>("face double", depth, new FaceVariable<NDIM, double>("face", depth));
    report_factory<double>("edge double", depth, new EdgeVariable<NDIM, double>("edge", depth));
    report_factory<double>("outerside double", depth, new OutersideVariable<NDIM, double>("outerside", depth));
    report_factory<double>("cell int queried as double", depth, new CellVariable<NDIM, int>("cell int", depth));
    report_factory<int>("cell int queried as int", depth, new CellVariable<NDIM, int>("cell int", depth));
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
    if (error == "factory")
    {
        Pointer<Variable<NDIM>> var = new CellVariable<NDIM, int>("unsupported factory");
        get_data_centering<double>(*var->getPatchDataFactory());
        return;
    }
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

using Cell = CartesianCentering<DataCentering::CELL>;
using Side = CartesianCentering<DataCentering::SIDE>;

VectorNd
position(const Patch<NDIM>& patch, const SAMRAI::hier::Index<NDIM>& index, const VectorNd& offset)
{
    const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch.getPatchGeometry();
    VectorNd x;
    for (int d = 0; d < NDIM; ++d)
    {
        x[d] = geometry->getXLower()[d] + geometry->getDx()[d] * (index(d) - patch.getBox().lower()(d) + offset[d]);
    }
    return x;
}

std::pair<double, double>
source_values(const VectorNd& x)
{
    double h = 0.25, liquid = 0.5;
    for (int d = 0; d < NDIM; ++d)
    {
        h += (d + 1) * x[d] * x[d] / 1024.0;
        liquid += (d + 2) * x[d] * x[d] / 2048.0;
    }
    return { h, liquid };
}

template <DataCentering C>
double
patchwise_error(PatchHierarchy<NDIM>& hierarchy, const int data_idx, const double time, const bool initial_time)
{
    using Layout = CartesianCentering<C>;
    double error = 0.0;
    for (int ln = 0; ln <= hierarchy.getFinestLevelNumber(); ++ln)
    {
        const Pointer<PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Pointer<typename Layout::template Data<double>> data = patch->getPatchData(data_idx);
            const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
            const double* const dx = geometry->getDx();
            for (int axis = 0; axis < Layout::num_axes(); ++axis)
            {
                const ArrayData<NDIM, double>& array = Layout::template array_data<double>(*data, axis);
                const Box<NDIM> interior = Layout::index_box(patch->getBox(), axis);
                double scale = 1.0;
                if constexpr (Layout::is_staggered())
                {
                    scale = dx[axis];
                }
                else
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        scale *= dx[d];
                    }
                }
                for (Box<NDIM>::Iterator it(array.getBox()); it; it++)
                {
                    const double actual = array(it(), 0);
                    if (!interior.contains(it()))
                    {
                        TBOX_ASSERT(std::isnan(actual));
                        continue;
                    }
                    const VectorNd x = position(*patch, it(), Layout::offset(axis));
                    std::pair<double, double> values = source_values(x);
                    if constexpr (Layout::is_staggered())
                    {
                        // The average of a quadratic at x +/- dx/2 includes this curvature term.
                        values.first += (axis + 1) * dx[axis] * dx[axis] / 4096.0;
                        values.second += (axis + 2) * dx[axis] * dx[axis] / 8192.0;
                    }
                    const double expected = 1.0 + 2.0 * values.first + 4.0 * values.first * values.second +
                                            time * scale + (initial_time ? 1.0 : 0.0);
                    TBOX_ASSERT(std::isfinite(actual));
                    error = std::max(error, std::abs(actual - expected));
                }
            }
        }
    }
    return IBTK_MPI::maxReduction(error);
}

// Compute a cell and a side field from cell sources with patchwise callbacks that capture state. The side
// field uses a mutable functor copied from a const lvalue, and the hierarchy and direct patch entry points are
// both exercised.
void
run_patchwise_cases(Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    VariableDatabase<NDIM>* const var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<VariableContext> context = var_db->getContext("patchwise");
    const Pointer<CellVariable<NDIM, double>> h_var = new CellVariable<NDIM, double>("H");
    const Pointer<CellVariable<NDIM, double>> liquid_var = new CellVariable<NDIM, double>("liquid fraction");
    const Pointer<CellVariable<NDIM, double>> cell_var = new CellVariable<NDIM, double>("cell blend");
    const Pointer<SideVariable<NDIM, double>> side_var = new SideVariable<NDIM, double>("side blend");
    IntVector<NDIM> output_ghosts(1);
    output_ghosts(1) = 2;
    const int h_idx = var_db->registerVariableAndContext(h_var, context, IntVector<NDIM>(1));
    const int liquid_idx = var_db->registerVariableAndContext(liquid_var, context, IntVector<NDIM>(1));
    const int cell_idx = var_db->registerVariableAndContext(cell_var, context, output_ghosts);
    const int side_idx = var_db->registerVariableAndContext(side_var, context, output_ghosts);
    TBOX_ASSERT(hierarchy->getFinestLevelNumber() == 1);

    int local_patches = 0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        const Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (const auto idx : { h_idx, liquid_idx, cell_idx, side_idx })
        {
            level->allocatePatchData(idx);
        }
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            ++local_patches;
            const Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Pointer<Cell::Data<double>> h = patch->getPatchData(h_idx);
            const Pointer<Cell::Data<double>> liquid = patch->getPatchData(liquid_idx);
            const Pointer<Cell::Data<double>> cell = patch->getPatchData(cell_idx);
            const Pointer<Side::Data<double>> side = patch->getPatchData(side_idx);
            cell->fillAll(std::numeric_limits<double>::quiet_NaN());
            side->fillAll(std::numeric_limits<double>::quiet_NaN());
            // Source ghosts are supplied by the caller, including physical and coarse-fine boundaries.
            for (Cell::Data<double>::Iterator it = Cell::begin(h->getGhostBox(), 0); it; it++)
            {
                const std::pair<double, double> values = source_values(position(*patch, it(), Cell::offset(0)));
                (*h)(it()) = values.first;
                (*liquid)(it()) = values.second;
            }
        }
    }
    TBOX_ASSERT(IBTK_MPI::sumReduction(local_patches) > 2);

    double expected_time = 0.25;
    bool expected_initial = true, expected_level = true;
    std::map<const Patch<NDIM>*, int> cell_calls, side_calls;
    const auto check_context =
        [&](Pointer<Patch<NDIM>> patch, const double time, const bool initial_time, Pointer<PatchLevel<NDIM>> level)
    {
        TBOX_ASSERT(time == expected_time && initial_time == expected_initial);
        TBOX_ASSERT(static_cast<bool>(level) == expected_level);
        if (level)
        {
            TBOX_ASSERT(level == hierarchy->getPatchLevel(patch->getPatchLevelNumber()));
            TBOX_ASSERT(level->getPatch(patch->getPatchNumber()) == patch);
        }
    };
    // A move-only functor, so the factory must take it by rvalue.
    auto cell_callback =
        [&, densities = std::make_unique<std::array<double, 3>>(std::array<double, 3>{ 1.0, 3.0, 7.0 })](
            const int data_idx,
            Pointer<Variable<NDIM>> var,
            Pointer<Patch<NDIM>> patch,
            const double time,
            const bool initial_time,
            Pointer<PatchLevel<NDIM>> level)
    {
        TBOX_ASSERT(data_idx == cell_idx && var == cell_var);
        check_context(patch, time, initial_time, level);
        ++cell_calls[patch.getPointer()];
        const Pointer<Cell::Data<double>> h_data = patch->getPatchData(h_idx);
        const Pointer<Cell::Data<double>> liquid_data = patch->getPatchData(liquid_idx);
        const Cell::Data<double>& h = *h_data;
        const Cell::Data<double>& liquid = *liquid_data;
        const Pointer<Cell::Data<double>> dst = patch->getPatchData(data_idx);
        const Pointer<CartesianPatchGeometry<NDIM>> patch_geometry = patch->getPatchGeometry();
        const double* const dx = patch_geometry->getDx();
        double volume = 1.0;
        for (int d = 0; d < NDIM; ++d)
        {
            volume *= dx[d];
        }
        for (Cell::Data<double>::Iterator it = Cell::begin(patch->getBox(), 0); it; it++)
        {
            const double material = (*densities)[1] + ((*densities)[2] - (*densities)[1]) * liquid(it());
            (*dst)(it()) =
                (*densities)[0] * (1.0 - h(it())) + material * h(it()) + time * volume + (initial_time ? 1.0 : 0.0);
        }
    };
    static_assert(PatchwiseCallback<decltype(cell_callback)>);
    static_assert(!std::copy_constructible<decltype(cell_callback)>);
    const Pointer<CartGridFunction> cell_function =
        make_cart_grid_patchwise_function("cell blend", std::move(cell_callback));

    Pointer<CartGridFunction> side_function;
    int observed_side_calls = 0;
    {
        // The factory must copy this const lvalue and invoke its mutable copy after this scope ends.
        const auto side_callback = [&, calls = 0](const int data_idx,
                                                  Pointer<Variable<NDIM>> var,
                                                  Pointer<Patch<NDIM>> patch,
                                                  const double time,
                                                  const bool initial_time,
                                                  Pointer<PatchLevel<NDIM>> level) mutable
        {
            TBOX_ASSERT(data_idx == side_idx && var == side_var);
            check_context(patch, time, initial_time, level);
            ++side_calls[patch.getPointer()];
            observed_side_calls = ++calls;
            const Pointer<Cell::Data<double>> h_data = patch->getPatchData(h_idx);
            const Pointer<Cell::Data<double>> liquid_data = patch->getPatchData(liquid_idx);
            const Cell::Data<double>& h = *h_data;
            const Cell::Data<double>& liquid = *liquid_data;
            const Pointer<Side::Data<double>> dst = patch->getPatchData(data_idx);
            const Pointer<CartesianPatchGeometry<NDIM>> patch_geometry = patch->getPatchGeometry();
            const double* const dx = patch_geometry->getDx();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Side::Data<double>::Iterator it = Side::begin(patch->getBox(), axis); it; it++)
                {
                    const Side::Index& side = it();
                    const double h_avg = 0.5 * (h(side.toCell(0)) + h(side.toCell(1)));
                    const double liquid_avg = 0.5 * (liquid(side.toCell(0)) + liquid(side.toCell(1)));
                    (*dst)(side) =
                        (1.0 - h_avg) + (3.0 + 4.0 * liquid_avg) * h_avg + time * dx[axis] + (initial_time ? 1.0 : 0.0);
                }
            }
        };
        static_assert(!PatchwiseCallback<decltype(side_callback)>);
        static_assert(PatchwiseCallback<std::decay_t<decltype(side_callback)>>);
        side_function = make_cart_grid_patchwise_function("side blend", side_callback);
    }
    TBOX_ASSERT(cell_function->isTimeDependent() && side_function->isTimeDependent());
    cell_function->setDataOnPatchHierarchy(cell_idx, cell_var, hierarchy, expected_time, expected_initial);
    const double cell_error =
        patchwise_error<DataCentering::CELL>(*hierarchy, cell_idx, expected_time, expected_initial);
    expected_time = 0.5;
    expected_initial = false;
    side_function->setDataOnPatchHierarchy(side_idx, side_var, hierarchy, expected_time, expected_initial);
    const double side_error =
        patchwise_error<DataCentering::SIDE>(*hierarchy, side_idx, expected_time, expected_initial);
    TBOX_ASSERT(observed_side_calls == local_patches);

    expected_time = 1.25;
    expected_level = false;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        const Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM>> patch = level->getPatch(p());
            TBOX_ASSERT(cell_calls[patch.getPointer()] == 1 && side_calls[patch.getPointer()] == 1);
            side_function->setDataOnPatch(side_idx, side_var, patch, expected_time);
            TBOX_ASSERT(side_calls[patch.getPointer()] == 2);
            const Pointer<Cell::Data<double>> h_data = patch->getPatchData(h_idx);
            const Pointer<Cell::Data<double>> liquid_data = patch->getPatchData(liquid_idx);
            const Cell::Data<double>& h = *h_data;
            const Cell::Data<double>& liquid = *liquid_data;
            for (Cell::Data<double>::Iterator it = Cell::begin(h.getGhostBox(), 0); it; it++)
            {
                const std::pair<double, double> expected = source_values(position(*patch, it(), Cell::offset(0)));
                TBOX_ASSERT(h(it()) == expected.first && liquid(it()) == expected.second);
            }
        }
    }
    TBOX_ASSERT(observed_side_calls == 2 * local_patches);
    TBOX_ASSERT(cell_calls.size() == static_cast<std::size_t>(local_patches));
    TBOX_ASSERT(side_calls.size() == static_cast<std::size_t>(local_patches));
    const double direct_error =
        patchwise_error<DataCentering::SIDE>(*hierarchy, side_idx, expected_time, expected_initial);
    plog << std::scientific << std::setprecision(12);
    plog << "cell blend error: " << cell_error << '\n';
    plog << "side blend errors (hierarchy, patch): " << side_error << ' ' << direct_error << '\n';
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
    if (app->getInputDatabase()->getBoolWithDefault("centering_helpers", false))
    {
        run_centering_helpers();
        return 0;
    }
    if (app->getInputDatabase()->getBoolWithDefault("patchwise", false))
    {
        run_patchwise_cases(hierarchy);
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
