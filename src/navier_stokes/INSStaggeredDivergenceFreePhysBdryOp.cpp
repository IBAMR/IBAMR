// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#include <ibamr/INSStaggeredDivergenceFreePhysBdryOp.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/StokesSpecifications.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/ShiftedPatchGeometry.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <BoxList.h>
#include <CartesianPatchGeometry.h>
#include <GridGeometry.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <Variable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "./ins_staggered_traction_stencil.h"

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
constexpr int REFINE_OP_STENCIL_WIDTH = 1;

// A value id names a value that the extension assigns or reads: id >= 0 is the temporary with that number, and id < 0
// is the face with index -id - 1 of the patch data.

// How a tangential velocity component is determined from the value across the boundary.
enum class TangentialRule
{
    PRESCRIBED,
    PSEUDO_TRACTION,
    TRACTION
};

// The location index of the boundary on the lower (upper) side of an axis.
inline unsigned int
location_index(const unsigned int axis, const bool is_upper)
{
    return 2 * axis + (is_upper ? 1 : 0);
}

// A set of boundaries: true at the location index of each boundary in the set.
using BoundarySet = std::array<bool, 2 * NDIM>;

inline int
count_boundaries(const BoundarySet& region)
{
    return static_cast<int>(std::count(region.begin(), region.end(), true));
}
} // namespace

/*!
 * A list of assignments in dependence order. Each assigns a weighted sum of values and boundary data to a temporary or
 * to a face of the patch data; a boundary datum is a value of g, divided by the viscosity for a traction condition.
 */
struct INSStaggeredDivergenceFreePhysBdryOp::Extension
{
    // A face of the patch data: the component and the offset in the data array.
    struct Face
    {
        unsigned int axis;
        int offset;
    };

    // weight times the value with the id value_id.
    struct Term
    {
        int value_id;
        double weight;
    };

    // weight times the value of a data point.
    struct DataTerm
    {
        int point;
        double weight;
    };

    // The assignment to the value id target of the terms and data terms that start at first_term and first_data_term.
    struct Assignment
    {
        int target;
        std::size_t first_term;
        std::size_t num_terms;
        std::size_t first_data_term;
        std::size_t num_data_terms;
    };

    // The value of the coefficient g of component axis at the boundary face index of the boundary location_index.
    struct DataPoint
    {
        unsigned int axis;
        unsigned int location_index;
        hier::Index<NDIM> index;
        bool divide_by_mu;
    };

    // The points evaluated by one call to the coefficient object.
    struct Batch
    {
        unsigned int axis;
        unsigned int location_index;
        Box<NDIM> box;
        std::vector<int> points;
    };

    // Compute the values of the data points on patch at fill_time from the coefficient g of bc_coefs; zero if
    // homogeneous_bc.
    std::vector<double> evaluateDataPoints(Patch<NDIM>& patch,
                                           const Pointer<Variable<NDIM>>& var,
                                           const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                           double fill_time,
                                           double mu,
                                           bool homogeneous_bc) const;

    // Assign the ghost values of the patch data u from its values in the domain and the values of the data points.
    void apply(const std::array<double*, NDIM>& u, const std::vector<double>& point_values) const;

    // Apply the transpose of apply() with zero data points to the patch data u.
    void applyTranspose(const std::array<double*, NDIM>& u) const;

    std::vector<Face> faces;
    int num_temporaries = 0;
    std::vector<Assignment> assignments;
    std::vector<Term> terms;
    std::vector<DataTerm> data_terms;
    std::vector<int> not_computable_faces;
    std::vector<DataPoint> points;
    std::vector<Batch> batches;
};

/*!
 * Construction of an Extension by the recursion that defines the ghost values.
 */
class INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder
{
public:
    // Prepare to build the extension of patch for the patch data data, filling the ghost values within
    // ghost_width_to_fill at fill_time. The patch must outlive the builder.
    ExtensionBuilder(const INSStaggeredHierarchyIntegrator& fluid_solver,
                     Patch<NDIM>& patch,
                     const Pointer<Variable<NDIM>>& var,
                     const SideData<NDIM, double>& data,
                     const IntVector<NDIM>& ghost_width_to_fill,
                     double fill_time);

    std::shared_ptr<Extension> build();

private:
    // A weighted sum of values and data points.
    struct LinearCombination
    {
        std::vector<Extension::Term> terms;
        std::vector<Extension::DataTerm> data_terms;
    };

    // The coefficients a and b of a component on a boundary.
    struct BcCoefficients
    {
        Box<NDIM> box;
        Pointer<ArrayData<NDIM, double>> acoef_data, bcoef_data;
    };

    // The key of a cached value: which value, the axis, and the index of the face.
    using Key = std::array<int, NDIM + 2>;

    static void add_term(std::vector<Extension::Term>& terms, int value_id, double weight);

    static void add_data_term(std::vector<Extension::DataTerm>& data_terms, int point, double weight);

    // Return the set of boundaries that the face i of component axis lies beyond.
    BoundarySet getRegion(unsigned int axis, const hier::Index<NDIM>& i) const;

    // Return the value id of the face i of component axis in the patch data, or nothing if it is not in the data.
    std::optional<int> getFaceValueId(unsigned int axis, const hier::Index<NDIM>& i);

    // Return the value id of the value of the face i of component axis, adding the assignments that compute it to the
    // extension, or nothing if it is not computable.
    std::optional<int> getValue(unsigned int axis, const hier::Index<NDIM>& i);

    // Return the value id of the value that the boundary loc gives to the face i, adding its assignments to the
    // extension, or nothing if it is not computable.
    std::optional<int> getValueAcross(unsigned int loc, unsigned int axis, const hier::Index<NDIM>& i);

    // Return the value id of the value of the face i of a cell in region: the value that the boundary loc gives to it
    // if it lies beyond the same boundaries as the cell, and the value it has otherwise.
    std::optional<int>
    getFieldValue(unsigned int loc, const BoundarySet& region, unsigned int axis, const hier::Index<NDIM>& i);

    // Set value to the value that the boundary loc gives to the face i of the component normal to it, which makes the
    // divergence of the adjacent ghost cell vanish; return false if it is not computable.
    bool setNormalValue(unsigned int loc, unsigned int axis, const hier::Index<NDIM>& i, LinearCombination& value);

    // Set value to the value that the boundary loc gives to the face i of a tangential component; return false if it is
    // not computable.
    bool setTangentialValue(unsigned int loc, unsigned int axis, const hier::Index<NDIM>& i, LinearCombination& value);

    // Add weight times the normal velocity at the face f of the boundary loc to value; return false if it is not
    // computable.
    bool addNormalVelocity(LinearCombination& value,
                           const hier::Index<NDIM>& f,
                           unsigned int loc,
                           unsigned int axis,
                           const BoundarySet& region,
                           double weight);

    // Return how component axis is determined across the boundary loc at the boundary face index.
    TangentialRule getRule(unsigned int axis, unsigned int loc, const hier::Index<NDIM>& index) const;

    // Return the index of the data point of component axis at the boundary face index of loc, adding it if needed.
    int getDataPoint(unsigned int axis, unsigned int loc, const hier::Index<NDIM>& index, bool divide_by_mu);

    // Evaluate the coefficients a and b of the boundary conditions.
    void readBcCoefficients();

    // Add the assignment of value to the value id target.
    void addAssignment(LinearCombination& value, int target);

    // Group the data points into batches that one call to a coefficient object evaluates.
    void groupDataPoints();

    int newTemporary();

    const INSStaggeredHierarchyIntegrator& d_fluid_solver;
    Patch<NDIM>& d_patch;
    Pointer<Variable<NDIM>> d_var;
    TractionBcType d_traction_bc_type;
    double d_fill_time;

    // The data and the geometry.
    std::array<Box<NDIM>, NDIM> d_data_box, d_fill_box;
    Box<NDIM> d_data_cell_box;
    std::array<double, NDIM> d_dx;
    std::array<int, NDIM> d_lower, d_upper;
    std::array<bool, NDIM> d_periodic;

    std::map<std::pair<unsigned int, unsigned int>, BcCoefficients> d_bc_coefficients;
    std::map<std::array<int, NDIM + 1>, int> d_face_ids;
    std::map<std::array<int, NDIM + 2>, int> d_point_ids;
    std::map<Key, std::optional<int>> d_cached_ids;
    std::shared_ptr<Extension> d_extension;
};

/////////////////////////////// EXTENSION BUILDER ////////////////////////////

INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::ExtensionBuilder(
    const INSStaggeredHierarchyIntegrator& fluid_solver,
    Patch<NDIM>& patch,
    const Pointer<Variable<NDIM>>& var,
    const SideData<NDIM, double>& data,
    const IntVector<NDIM>& ghost_width_to_fill,
    const double fill_time)
    : d_fluid_solver(fluid_solver),
      d_patch(patch),
      d_var(var),
      d_traction_bc_type(fluid_solver.getTractionBcType()),
      d_fill_time(fill_time),
      d_extension(std::make_shared<Extension>())
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(data.getGhostCellWidth(), ghost_width_to_fill);
    d_data_cell_box = data.getGhostBox();
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        d_data_box[axis] = data.getArrayData(axis).getBox();
        d_fill_box[axis] = SideGeometry<NDIM>::toSideBox(Box<NDIM>::grow(patch.getBox(), gcw_to_fill), axis);
        d_dx[axis] = pgeom->getDx()[axis];
    }

    // The domain must be a single box.
    const Pointer<PatchHierarchy<NDIM>> hierarchy = fluid_solver.getPatchHierarchy();
    const BoxArray<NDIM> domain = traction_stencil::get_physical_domain(hierarchy, pgeom->getRatio());
    Box<NDIM> domain_box = domain[0];
    for (int n = 1; n < domain.size(); ++n)
    {
        domain_box += domain[n];
    }
    BoxList<NDIM> uncovered(domain_box);
    uncovered.removeIntersections(BoxList<NDIM>(domain));
    if (uncovered.getNumberOfItems() != 0)
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::ExtensionBuilder():\n"
                   << "  the physical domain must be a single rectangular box" << std::endl);
    }
    const IntVector<NDIM> periodic_shift =
        fluid_solver.getPatchHierarchy()->getGridGeometry()->getPeriodicShift(pgeom->getRatio());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_lower[d] = domain_box.lower(d);
        d_upper[d] = domain_box.upper(d);
        d_periodic[d] = periodic_shift(d) != 0;
        if (!d_periodic[d] && d_upper[d] - d_lower[d] + 1 < data.getGhostCellWidth().max())
        {
            TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::ExtensionBuilder():\n"
                       << "  the physical domain is " << d_upper[d] - d_lower[d] + 1 << " cells wide along axis " << d
                       << ", but the ghost cell width of the velocity is " << data.getGhostCellWidth().max() << "\n"
                       << "  the domain must be at least as wide as the ghost cell width" << std::endl);
        }
    }
}

std::shared_ptr<INSStaggeredDivergenceFreePhysBdryOp::Extension>
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::build()
{
    readBcCoefficients();
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (Box<NDIM>::Iterator b(d_fill_box[axis]); b; b++)
        {
            const hier::Index<NDIM>& i = b();
            if (count_boundaries(getRegion(axis, i)) != 0 && !getValue(axis, i))
            {
                d_extension->not_computable_faces.push_back(-*getFaceValueId(axis, i) - 1);
            }
        }
    }
    groupDataPoints();
    return d_extension;
}

BoundarySet
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::getRegion(const unsigned int axis,
                                                                  const hier::Index<NDIM>& i) const
{
    BoundarySet region{};
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d_periodic[d])
        {
            continue;
        }
        if (i(d) < d_lower[d])
        {
            region[location_index(d, false)] = true;
        }
        else if (i(d) > d_upper[d] + (d == axis ? 1 : 0))
        {
            region[location_index(d, true)] = true;
        }
    }
    return region;
}

std::optional<int>
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::getFaceValueId(const unsigned int axis,
                                                                       const hier::Index<NDIM>& i)
{
    if (!d_data_box[axis].contains(i))
    {
        return std::nullopt;
    }
    std::array<int, NDIM + 1> key;
    key[0] = static_cast<int>(axis);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        key[d + 1] = i(d);
    }
    const auto it = d_face_ids.find(key);
    if (it != d_face_ids.end())
    {
        return -it->second - 1;
    }
    const int id = static_cast<int>(d_extension->faces.size());
    d_extension->faces.push_back({ axis, d_data_box[axis].offset(i) });
    d_face_ids[key] = id;
    return -id - 1;
}

int
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::newTemporary()
{
    return d_extension->num_temporaries++;
}

void
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::add_term(std::vector<Extension::Term>& terms,
                                                                 const int value_id,
                                                                 const double weight)
{
    for (auto& term : terms)
    {
        if (term.value_id == value_id)
        {
            term.weight += weight;
            return;
        }
    }
    terms.push_back({ value_id, weight });
}

void
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::add_data_term(std::vector<Extension::DataTerm>& data_terms,
                                                                      const int point,
                                                                      const double weight)
{
    for (auto& term : data_terms)
    {
        if (term.point == point)
        {
            term.weight += weight;
            return;
        }
    }
    data_terms.push_back({ point, weight });
}

std::optional<int>
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::getValue(const unsigned int axis, const hier::Index<NDIM>& i)
{
    if (!d_data_box[axis].contains(i))
    {
        return std::nullopt;
    }
    const BoundarySet region = getRegion(axis, i);
    if (count_boundaries(region) == 0)
    {
        return getFaceValueId(axis, i);
    }
    Key key;
    key[0] = 0;
    key[1] = static_cast<int>(axis);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        key[d + 2] = i(d);
    }
    const auto it = d_cached_ids.find(key);
    if (it != d_cached_ids.end())
    {
        return it->second;
    }

    // Faces that are stored are assigned in the patch data; the others get temporaries.
    const int target = d_fill_box[axis].contains(i) ? *getFaceValueId(axis, i) : newTemporary();
    const int n = count_boundaries(region);
    LinearCombination value;
    bool computable = true;
    for (unsigned int loc = 0; loc < 2 * NDIM && computable; ++loc)
    {
        if (!region[loc])
        {
            continue;
        }
        if (n == 1)
        {
            computable =
                (loc / 2 == axis) ? setNormalValue(loc, axis, i, value) : setTangentialValue(loc, axis, i, value);
        }
        else
        {
            const std::optional<int> across_id = getValueAcross(loc, axis, i);
            computable = across_id.has_value();
            if (computable)
            {
                add_term(value.terms, *across_id, 1.0 / n);
            }
        }
    }
    std::optional<int> result;
    if (computable)
    {
        addAssignment(value, target);
        result = target;
    }
    d_cached_ids[key] = result;
    return result;
}

std::optional<int>
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::getValueAcross(const unsigned int loc,
                                                                       const unsigned int axis,
                                                                       const hier::Index<NDIM>& i)
{
    const BoundarySet region = getRegion(axis, i);
    if (count_boundaries(region) == 1)
    {
        return getValue(axis, i);
    }
    if (!d_data_box[axis].contains(i))
    {
        return std::nullopt;
    }
    Key key;
    key[0] = 1 + static_cast<int>(loc);
    key[1] = static_cast<int>(axis);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        key[d + 2] = i(d);
    }
    const auto it = d_cached_ids.find(key);
    if (it != d_cached_ids.end())
    {
        return it->second;
    }
    const int target = newTemporary();
    LinearCombination value;
    const bool computable =
        (loc / 2 == axis) ? setNormalValue(loc, axis, i, value) : setTangentialValue(loc, axis, i, value);
    std::optional<int> result;
    if (computable)
    {
        addAssignment(value, target);
        result = target;
    }
    d_cached_ids[key] = result;
    return result;
}

std::optional<int>
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::getFieldValue(const unsigned int loc,
                                                                      const BoundarySet& region,
                                                                      const unsigned int axis,
                                                                      const hier::Index<NDIM>& i)
{
    return getRegion(axis, i) == region ? getValueAcross(loc, axis, i) : getValue(axis, i);
}

bool
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::setNormalValue(const unsigned int loc,
                                                                       const unsigned int axis,
                                                                       const hier::Index<NDIM>& i,
                                                                       LinearCombination& value)
{
    const bool is_upper = loc % 2 == 1;
    const BoundarySet region = getRegion(axis, i);
    hier::Index<NDIM> cell = i;
    hier::Index<NDIM> inner = i;
    if (is_upper)
    {
        cell(axis) -= 1;
        inner(axis) -= 1;
    }
    else
    {
        inner(axis) += 1;
    }
    const std::optional<int> inner_id = getFieldValue(loc, region, axis, inner);
    if (!inner_id)
    {
        return false;
    }
    add_term(value.terms, *inner_id, 1.0);
    for (unsigned int t = 0; t < NDIM; ++t)
    {
        if (t == axis)
        {
            continue;
        }
        hier::Index<NDIM> cell_hi = cell;
        cell_hi(t) += 1;
        const std::optional<int> hi_id = getFieldValue(loc, region, t, cell_hi);
        const std::optional<int> lo_id = getFieldValue(loc, region, t, cell);
        if (!hi_id || !lo_id)
        {
            return false;
        }
        const double weight = (is_upper ? -1.0 : 1.0) * d_dx[axis] / d_dx[t];
        add_term(value.terms, *hi_id, weight);
        add_term(value.terms, *lo_id, -weight);
    }
    return true;
}

bool
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::setTangentialValue(const unsigned int loc,
                                                                           const unsigned int axis,
                                                                           const hier::Index<NDIM>& i,
                                                                           LinearCombination& value)
{
    const unsigned int normal_axis = loc / 2;
    const bool is_upper = loc % 2 == 1;
    const BoundarySet region = getRegion(axis, i);
    const int k = is_upper ? i(normal_axis) - d_upper[normal_axis] - 1 : d_lower[normal_axis] - 1 - i(normal_axis);
    hier::Index<NDIM> mirror = i;
    mirror(normal_axis) = is_upper ? d_upper[normal_axis] - k : d_lower[normal_axis] + k;
    const std::optional<int> mirror_id = getValue(axis, mirror);
    if (!mirror_id)
    {
        return false;
    }
    hier::Index<NDIM> boundary_face = i;
    boundary_face(normal_axis) = is_upper ? d_upper[normal_axis] + 1 : d_lower[normal_axis];
    const double sgn = is_upper ? +1.0 : -1.0;
    const double gamma_weight = (2.0 * k + 1.0) * d_dx[normal_axis] * sgn;
    switch (getRule(axis, loc, boundary_face))
    {
    case TangentialRule::PRESCRIBED:
        add_term(value.terms, *mirror_id, -1.0);
        add_data_term(value.data_terms, getDataPoint(axis, loc, boundary_face, false), 2.0);
        break;
    case TangentialRule::PSEUDO_TRACTION:
        add_term(value.terms, *mirror_id, 1.0);
        add_data_term(value.data_terms, getDataPoint(axis, loc, boundary_face, true), gamma_weight);
        break;
    case TangentialRule::TRACTION:
    {
        add_term(value.terms, *mirror_id, 1.0);
        add_data_term(value.data_terms, getDataPoint(axis, loc, boundary_face, true), gamma_weight);
        hier::Index<NDIM> boundary_face_below = boundary_face;
        boundary_face_below(axis) -= 1;
        const double weight = gamma_weight / d_dx[axis];
        if (!addNormalVelocity(value, boundary_face, loc, axis, region, -weight) ||
            !addNormalVelocity(value, boundary_face_below, loc, axis, region, +weight))
        {
            return false;
        }
        break;
    }
    }
    return true;
}

bool
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::addNormalVelocity(LinearCombination& value,
                                                                          const hier::Index<NDIM>& f,
                                                                          const unsigned int loc,
                                                                          const unsigned int axis,
                                                                          const BoundarySet& region,
                                                                          const double weight)
{
    const unsigned int normal_axis = loc / 2;
    BoundarySet other_region = region;
    other_region[loc] = false;
    const BoundarySet f_region = getRegion(normal_axis, f);
    BoundarySet extra{};
    for (unsigned int l = 0; l < 2 * NDIM; ++l)
    {
        extra[l] = f_region[l] && !other_region[l];
    }

    // A face that is not beyond a corner has its value. A face one cell beyond a corner, across a boundary not in the
    // region, has the corner value: 2*u_b - u if the adjacent boundary prescribes the normal velocity u_b, with u the
    // nearest face on B, and otherwise the linear extrapolation along B.
    if (count_boundaries(extra) == 0)
    {
        const std::optional<int> id = getValue(normal_axis, f);
        if (!id)
        {
            return false;
        }
        add_term(value.terms, *id, weight);
        return true;
    }
    const unsigned int adj_loc_lower = location_index(axis, false);
    const unsigned int adj_loc_upper = location_index(axis, true);
    const auto is_only_extra = [&](const unsigned int l) { return extra[l] && count_boundaries(extra) == 1; };
    if (!is_only_extra(adj_loc_lower) && !is_only_extra(adj_loc_upper))
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::addNormalVelocity():\n"
                   << "  the normal velocity face " << f << " is not one cell beyond a corner of the boundary " << loc
                   << std::endl);
    }
    const int step = is_only_extra(adj_loc_lower) ? +1 : -1;

    // A boundary segment of one face has no second face to extrapolate from.
    const int j_next = f(axis) + 2 * step;
    const bool can_extrapolate = j_next >= d_lower[axis] && j_next <= d_upper[axis];
    traction_stencil::NormalVelocityStencil stencil =
        traction_stencil::get_corner_normal_velocity_stencil(f, normal_axis, axis, step, can_extrapolate);
    hier::Index<NDIM> adjacent_boundary_face = f;
    adjacent_boundary_face(axis) = stencil.adjacent_location_index % 2 == 1 ? d_upper[axis] + 1 : d_lower[axis];
    if (getRule(normal_axis, stencil.adjacent_location_index, adjacent_boundary_face) == TangentialRule::PRESCRIBED)
    {
        stencil.useAdjacentBoundaryValue();
    }
    for (int n = 0; n < 2; ++n)
    {
        if (stencil.weight[n] == 0.0)
        {
            continue;
        }
        const std::optional<int> id = getValue(normal_axis, stencil.idx[n]);
        if (!id)
        {
            return false;
        }
        add_term(value.terms, *id, weight * stencil.weight[n]);
    }
    if (stencil.uses_adjacent_boundary_value)
    {
        add_data_term(value.data_terms,
                      getDataPoint(normal_axis, stencil.adjacent_location_index, adjacent_boundary_face, false),
                      2.0 * weight);
    }
    return true;
}

void
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::readBcCoefficients()
{
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs = d_fluid_solver.getPhysicalBoundaryConditions();
    Box<NDIM> cell_box = d_data_cell_box;
    cell_box.grow(1);

    // Evaluate the coefficients of a component at its faces, with the patch geometry shifted by half a cell.
    std::array<std::optional<ShiftedPatchGeometry>, NDIM> shifted_geometry;
    for (unsigned int normal_axis = 0; normal_axis < NDIM; ++normal_axis)
    {
        if (d_periodic[normal_axis])
        {
            continue;
        }
        for (int side = 0; side < 2; ++side)
        {
            const bool is_upper = side == 1;
            const unsigned int loc = location_index(normal_axis, is_upper);
            if (is_upper ? cell_box.upper(normal_axis) <= d_upper[normal_axis] :
                           cell_box.lower(normal_axis) >= d_lower[normal_axis])
            {
                continue;
            }
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                if (axis == normal_axis)
                {
                    continue;
                }
                if (bc_coefs.size() != NDIM || !bc_coefs[axis])
                {
                    TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::readBcCoefficients():\n"
                               << "  the fluid solver has no velocity boundary conditions for component " << axis
                               << "\n"
                               << "  register them with registerPhysicalBoundaryConditions()" << std::endl);
                }
                BcCoefficients bc_coefficients;
                bc_coefficients.box = SideGeometry<NDIM>::toSideBox(cell_box, axis);
                const int plane = is_upper ? d_upper[normal_axis] + 1 : d_lower[normal_axis];
                bc_coefficients.box.lower()(normal_axis) = plane;
                bc_coefficients.box.upper()(normal_axis) = plane;
                bc_coefficients.acoef_data = new ArrayData<NDIM, double>(bc_coefficients.box, 1);
                bc_coefficients.bcoef_data = new ArrayData<NDIM, double>(bc_coefficients.box, 1);
                Pointer<ArrayData<NDIM, double>> gcoef_data;
                if (!shifted_geometry[axis])
                {
                    std::array<double, NDIM> shift;
                    shift.fill(0.0);
                    shift[axis] = -0.5 * d_dx[axis];
                    shifted_geometry[axis].emplace(d_patch, shift);
                }
                Box<NDIM> bdry_cell_box = bc_coefficients.box;
                if (!is_upper)
                {
                    bdry_cell_box.shift(normal_axis, -1);
                }
                const BoundaryBox<NDIM> bdry_box(bdry_cell_box, 1, loc);
                {
                    const ShiftedPatchGeometry::Scope scope(d_patch, *shifted_geometry[axis]);
                    bc_coefs[axis]->setBcCoefs(bc_coefficients.acoef_data,
                                               bc_coefficients.bcoef_data,
                                               gcoef_data,
                                               d_var,
                                               d_patch,
                                               bdry_box,
                                               d_fill_time);
                }
                d_bc_coefficients[std::make_pair(axis, loc)] = bc_coefficients;
            }
        }
    }
}

TangentialRule
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::getRule(const unsigned int axis,
                                                                const unsigned int loc,
                                                                const hier::Index<NDIM>& index) const
{
    const auto it = d_bc_coefficients.find(std::make_pair(axis, loc));
    if (it == d_bc_coefficients.end() || !it->second.box.contains(index))
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::getRule():\n"
                   << "  the boundary condition of component " << axis << " on the boundary " << loc
                   << " is not known at the boundary face " << index << std::endl);
    }
    const double alpha = (*it->second.acoef_data)(index, 0);
    const double beta = (*it->second.bcoef_data)(index, 0);
    const bool velocity_bc = (alpha == 1.0 && beta == 0.0);
    const bool traction_bc = (alpha == 0.0 && beta == 1.0);
    if (!velocity_bc && !traction_bc)
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::getRule():\n"
                   << "  the boundary condition of component " << axis << " on the boundary " << loc
                   << " at the boundary face " << index << " has unsupported coefficients (a, b) = (" << alpha << ", "
                   << beta << ")\n"
                   << "  use (a, b) = (1, 0) to prescribe the velocity or (0, 1) to prescribe the traction"
                   << std::endl);
    }
    if (velocity_bc)
    {
        return TangentialRule::PRESCRIBED;
    }
    switch (d_traction_bc_type)
    {
    case TRACTION:
        return TangentialRule::TRACTION;
    case PSEUDO_TRACTION:
        return TangentialRule::PSEUDO_TRACTION;
    default:
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::getRule():\n"
                   << "  unrecognized or unsupported traction boundary condition type: "
                   << enum_to_string<TractionBcType>(d_traction_bc_type) << "\n"
                   << "  use TRACTION or PSEUDO_TRACTION" << std::endl);
    }
    return TangentialRule::TRACTION;
}

int
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::getDataPoint(const unsigned int axis,
                                                                     const unsigned int loc,
                                                                     const hier::Index<NDIM>& index,
                                                                     const bool divide_by_mu)
{
    std::array<int, NDIM + 2> key;
    key[0] = static_cast<int>(axis);
    key[1] = static_cast<int>(loc);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        key[d + 2] = index(d);
    }
    const auto it = d_point_ids.find(key);
    if (it != d_point_ids.end())
    {
        return it->second;
    }
    const int id = static_cast<int>(d_extension->points.size());
    d_extension->points.push_back({ axis, loc, index, divide_by_mu });
    d_point_ids[key] = id;
    return id;
}

void
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::addAssignment(LinearCombination& value, const int target)
{
    Extension::Assignment assignment;
    assignment.target = target;
    assignment.first_term = d_extension->terms.size();
    assignment.num_terms = value.terms.size();
    assignment.first_data_term = d_extension->data_terms.size();
    assignment.num_data_terms = value.data_terms.size();
    d_extension->terms.insert(d_extension->terms.end(), value.terms.begin(), value.terms.end());
    d_extension->data_terms.insert(d_extension->data_terms.end(), value.data_terms.begin(), value.data_terms.end());
    d_extension->assignments.push_back(assignment);
}

void
INSStaggeredDivergenceFreePhysBdryOp::ExtensionBuilder::groupDataPoints()
{
    std::map<std::pair<unsigned int, unsigned int>, std::size_t> batch_ids;
    for (int p = 0; p < static_cast<int>(d_extension->points.size()); ++p)
    {
        const Extension::DataPoint& point = d_extension->points[p];
        const auto key = std::make_pair(point.axis, point.location_index);
        auto it = batch_ids.find(key);
        if (it == batch_ids.end())
        {
            it = batch_ids.insert(std::make_pair(key, d_extension->batches.size())).first;
            d_extension->batches.push_back(
                { point.axis, point.location_index, Box<NDIM>(point.index, point.index), {} });
        }
        Extension::Batch& batch = d_extension->batches[it->second];
        batch.box += Box<NDIM>(point.index, point.index);
        batch.points.push_back(p);
    }
}

/////////////////////////////// TAPE /////////////////////////////////////////

std::vector<double>
INSStaggeredDivergenceFreePhysBdryOp::Extension::evaluateDataPoints(
    Patch<NDIM>& patch,
    const Pointer<Variable<NDIM>>& var,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    const double fill_time,
    const double mu,
    const bool homogeneous_bc) const
{
    std::vector<double> values(points.size(), 0.0);
    if (homogeneous_bc)
    {
        return values;
    }
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    std::array<std::optional<ShiftedPatchGeometry>, NDIM> shifted_geometry;
    for (const auto& batch : batches)
    {
        if (!shifted_geometry[batch.axis])
        {
            std::array<double, NDIM> shift;
            shift.fill(0.0);
            shift[batch.axis] = -0.5 * pgeom->getDx()[batch.axis];
            shifted_geometry[batch.axis].emplace(patch, shift);
        }
        Box<NDIM> bdry_cell_box = batch.box;
        if (batch.location_index % 2 == 0)
        {
            bdry_cell_box.shift(batch.location_index / 2, -1);
        }
        const BoundaryBox<NDIM> bdry_box(bdry_cell_box, 1, batch.location_index);
        Pointer<ArrayData<NDIM, double>> acoef_data, bcoef_data;
        Pointer<ArrayData<NDIM, double>> gcoef_data = new ArrayData<NDIM, double>(batch.box, 1);
        {
            const ShiftedPatchGeometry::Scope scope(patch, *shifted_geometry[batch.axis]);
            bc_coefs[batch.axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, patch, bdry_box, fill_time);
        }
        for (const int p : batch.points)
        {
            const double g = (*gcoef_data)(points[p].index, 0);
            values[p] = points[p].divide_by_mu ? g / mu : g;
        }
    }
    return values;
}

void
INSStaggeredDivergenceFreePhysBdryOp::Extension::apply(const std::array<double*, NDIM>& u,
                                                       const std::vector<double>& point_values) const
{
    std::vector<double> temporaries(num_temporaries, 0.0);
    const auto value = [&](const int value_id) -> double&
    {
        if (value_id >= 0)
        {
            return temporaries[value_id];
        }
        const Face& face = faces[-value_id - 1];
        return u[face.axis][face.offset];
    };
    for (const auto& assignment : assignments)
    {
        double sum = 0.0;
        for (std::size_t n = assignment.first_term; n < assignment.first_term + assignment.num_terms; ++n)
        {
            sum += terms[n].weight * value(terms[n].value_id);
        }
        for (std::size_t n = assignment.first_data_term; n < assignment.first_data_term + assignment.num_data_terms;
             ++n)
        {
            sum += data_terms[n].weight * point_values[data_terms[n].point];
        }
        value(assignment.target) = sum;
    }
    for (const int id : not_computable_faces)
    {
        u[faces[id].axis][faces[id].offset] = std::numeric_limits<double>::quiet_NaN();
    }
}

void
INSStaggeredDivergenceFreePhysBdryOp::Extension::applyTranspose(const std::array<double*, NDIM>& u) const
{
    std::vector<double> temporaries(num_temporaries, 0.0);
    const auto value = [&](const int value_id) -> double&
    {
        if (value_id >= 0)
        {
            return temporaries[value_id];
        }
        const Face& face = faces[-value_id - 1];
        return u[face.axis][face.offset];
    };
    for (auto it = assignments.rbegin(); it != assignments.rend(); ++it)
    {
        const double accumulated = value(it->target);
        value(it->target) = 0.0;
        for (std::size_t n = it->first_term; n < it->first_term + it->num_terms; ++n)
        {
            value(terms[n].value_id) += terms[n].weight * accumulated;
        }
    }
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredDivergenceFreePhysBdryOp::INSStaggeredDivergenceFreePhysBdryOp(
    const INSStaggeredHierarchyIntegrator* const fluid_solver,
    const bool homogeneous_bc)
    : d_fluid_solver(fluid_solver)
{
    setHomogeneousBc(homogeneous_bc);
    return;
} // INSStaggeredDivergenceFreePhysBdryOp

INSStaggeredDivergenceFreePhysBdryOp::~INSStaggeredDivergenceFreePhysBdryOp() = default;

void
INSStaggeredDivergenceFreePhysBdryOp::clearCache()
{
    d_extension_cache.clear();
    return;
} // clearCache

void
INSStaggeredDivergenceFreePhysBdryOp::setPhysicalBcCoef(RobinBcCoefStrategy<NDIM>* /*bc_coef*/)
{
    TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::setPhysicalBcCoef():\n"
               << "  this class uses the velocity boundary conditions of the fluid solver\n"
               << "  register them with registerPhysicalBoundaryConditions() of the fluid solver" << std::endl);
    return;
} // setPhysicalBcCoef

void
INSStaggeredDivergenceFreePhysBdryOp::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& /*bc_coefs*/)
{
    TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::setPhysicalBcCoefs():\n"
               << "  this class uses the velocity boundary conditions of the fluid solver\n"
               << "  register them with registerPhysicalBoundaryConditions() of the fluid solver" << std::endl);
    return;
} // setPhysicalBcCoefs

void
INSStaggeredDivergenceFreePhysBdryOp::setPhysicalBoundaryConditions(Patch<NDIM>& patch,
                                                                    const double fill_time,
                                                                    const IntVector<NDIM>& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector<NDIM>(0))
    {
        return;
    }
    const double mu = d_fluid_solver->getStokesSpecifications()->getMu();
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        const std::shared_ptr<const Extension> extension =
            getExtension(patch, patch_data_idx, ghost_width_to_fill, fill_time);
        Pointer<SideData<NDIM, double>> data = patch.getPatchData(patch_data_idx);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<Variable<NDIM>> var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        const std::vector<double> point_values = extension->evaluateDataPoints(
            patch, var, d_fluid_solver->getPhysicalBoundaryConditions(), fill_time, mu, d_homogeneous_bc);
        std::array<double*, NDIM> u;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            u[axis] = data->getPointer(axis, 0);
        }
        extension->apply(u, point_values);
    }
    return;
} // setPhysicalBoundaryConditions

IntVector<NDIM>
INSStaggeredDivergenceFreePhysBdryOp::getRefineOpStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
INSStaggeredDivergenceFreePhysBdryOp::accumulateFromPhysicalBoundaryData(Patch<NDIM>& patch,
                                                                         const double fill_time,
                                                                         const IntVector<NDIM>& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector<NDIM>(0))
    {
        return;
    }
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        const std::shared_ptr<const Extension> extension =
            getExtension(patch, patch_data_idx, ghost_width_to_fill, fill_time);
        Pointer<SideData<NDIM, double>> data = patch.getPatchData(patch_data_idx);
        std::array<double*, NDIM> u;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            u[axis] = data->getPointer(axis, 0);
        }
        extension->applyTranspose(u);
    }
    return;
} // accumulateFromPhysicalBoundaryData

/////////////////////////////// PRIVATE //////////////////////////////////////

std::shared_ptr<const INSStaggeredDivergenceFreePhysBdryOp::Extension>
INSStaggeredDivergenceFreePhysBdryOp::getExtension(Patch<NDIM>& patch,
                                                   const int patch_data_idx,
                                                   const IntVector<NDIM>& ghost_width_to_fill,
                                                   const double fill_time)
{
    Pointer<SideData<NDIM, double>> data = patch.getPatchData(patch_data_idx);
    if (!data)
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::getExtension():\n"
                   << "  patch data index " << patch_data_idx
                   << " does not correspond to a side-centered double precision variable" << std::endl);
    }
    if (data->getDepth() != 1)
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::getExtension():\n"
                   << "  the depth of patch data index " << patch_data_idx << " is " << data->getDepth()
                   << ", but a velocity has depth 1" << std::endl);
    }
    const IntVector<NDIM> gcw_to_fill = IntVector<NDIM>::min(data->getGhostCellWidth(), ghost_width_to_fill);

    // Patches of the hierarchy are identified by their level and box.
    Pointer<PatchHierarchy<NDIM>> hierarchy = d_fluid_solver->getPatchHierarchy();
    const int ln = patch.getPatchLevelNumber();
    bool in_hierarchy = false;
    if (hierarchy && ln >= 0 && ln <= hierarchy->getFinestLevelNumber())
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        const int pn = patch.getPatchNumber();
        in_hierarchy = pn >= 0 && pn < level->getNumberOfPatches() && level->getPatch(pn).getPointer() == &patch;
    }
    std::array<int, 4 * NDIM + 1> key;
    key[0] = ln;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        key[1 + d] = patch.getBox().lower(d);
        key[1 + NDIM + d] = patch.getBox().upper(d);
        key[1 + 2 * NDIM + d] = data->getGhostCellWidth()(d);
        key[1 + 3 * NDIM + d] = gcw_to_fill(d);
    }
    if (in_hierarchy)
    {
        const auto it = d_extension_cache.find(key);
        if (it != d_extension_cache.end())
        {
            return it->second;
        }
    }
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM>> var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    ExtensionBuilder builder(*d_fluid_solver, patch, var, *data, gcw_to_fill, fill_time);
    std::shared_ptr<const Extension> extension = builder.build();
    if (in_hierarchy)
    {
        d_extension_cache[key] = extension;
    }
    return extension;
} // getExtension

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
