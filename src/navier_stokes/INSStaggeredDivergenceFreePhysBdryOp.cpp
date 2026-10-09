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

#include <ibtk/ibtk_utilities.h>

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

#include <array>
#include <limits>
#include <map>
#include <memory>
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

// A slot is a value that the tape assigns or reads. A slot s >= 0 is the
// temporary value with index s. A slot s < 0 is the face with index -s - 1 of
// the patch data. A value that is not computable has the slot NOT_COMPUTABLE.
constexpr int NOT_COMPUTABLE = std::numeric_limits<int>::min();

// How a velocity component tangential to a boundary is determined from the
// value across the boundary.
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

// The number of boundaries in the set region.
inline int
count_boundaries(const unsigned int region)
{
    int n = 0;
    for (unsigned int r = region; r != 0; r >>= 1)
    {
        n += static_cast<int>(r & 1u);
    }
    return n;
}
} // namespace

/*!
 * The assignments of a Tape are in an order in which each value follows the
 * values it depends on. The assignment to a slot is the sum over terms of
 * weight * value of slot, plus the sum over data terms of weight * value of data
 * point. Data points are values of the coefficient g, divided by the viscosity
 * for a traction condition.
 */
struct INSStaggeredDivergenceFreePhysBdryOp::Tape
{
    // A face of the patch data: the component and the offset in the data array.
    struct Face
    {
        unsigned int axis;
        int offset;
    };

    // weight times the value of a slot.
    struct Term
    {
        int slot;
        double weight;
    };

    // weight times the value of a data point.
    struct DataTerm
    {
        int point;
        double weight;
    };

    // The assignment to the slot target of the sum of the terms and data terms
    // that start at first_term and first_data_term.
    struct Entry
    {
        int target;
        std::size_t first_term;
        std::size_t num_terms;
        std::size_t first_data_term;
        std::size_t num_data_terms;
    };

    // The value of the coefficient g of component axis at the boundary face
    // with index index on the boundary with the given location index.
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

    /*!
     * Compute the values of the data points of the tape on patch at
     * fill_time by evaluating the coefficient g of the boundary conditions
     * bc_coefs. The coefficients of a tangential component are evaluated at its
     * faces on the boundary. The values are zero if homogeneous_bc is true.
     */
    std::vector<double> evaluateDataPoints(Patch<NDIM>& patch,
                                           const Pointer<Variable<NDIM>>& var,
                                           const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                           double fill_time,
                                           double mu,
                                           bool homogeneous_bc) const;

    /*!
     * Assign the ghost values in the patch data with the pointers u to the
     * values in the patch data and the values of the data points.
     */
    void apply(const std::array<double*, NDIM>& u, const std::vector<double>& point_values) const;

    /*!
     * Apply the transpose of apply() with zero data points to the values in
     * the patch data with the pointers u.
     */
    void applyTranspose(const std::array<double*, NDIM>& u) const;

    /*!
     * Return the index in faces of a face whose ghost value is not computable and whose value in the patch data with
     * the pointers u is not zero, or -1 if there is none.
     */
    int findNonzeroNotComputableFace(const std::array<double*, NDIM>& u) const;

    std::vector<Face> faces;
    int num_temporaries = 0;
    std::vector<Entry> entries;
    std::vector<Term> terms;
    std::vector<DataTerm> data_terms;
    std::vector<int> not_computable_faces;
    std::vector<DataPoint> points;
    std::vector<Batch> batches;
};

/*!
 * Computation of a Tape by the recursion that defines the ghost values.
 */
class INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder
{
public:
    /*!
     * Prepare to build the tape of patch for the patch data data with the ghost
     * box of data, filling the ghost values within ghost_width_to_fill at
     * fill_time. The patch must outlive the builder. Its geometry is replaced
     * temporarily while the boundary coefficients are evaluated.
     */
    TapeBuilder(const INSStaggeredHierarchyIntegrator& fluid_solver,
                traction_stencil::PhysicalDomainCache& domain_cache,
                Patch<NDIM>& patch,
                const Pointer<Variable<NDIM>>& var,
                const SideData<NDIM, double>& data,
                const IntVector<NDIM>& ghost_width_to_fill,
                double fill_time);

    /*!
     * Compute and return the tape.
     */
    std::shared_ptr<Tape> build();

private:
    // A weighted sum of slots and data points.
    struct Combination
    {
        std::vector<Tape::Term> terms;
        std::vector<Tape::DataTerm> data_terms;
    };

    // The coefficients a and b of a component on a boundary.
    struct Rules
    {
        Box<NDIM> box;
        Pointer<ArrayData<NDIM, double>> acoef_data, bcoef_data;
    };

    // The key of a memoized value: which value is computed, the axis, and the
    // index of the face.
    using Key = std::array<int, NDIM + 2>;

    /*!
     * Add weight times slot to the sum of terms.
     */
    static void add_term(std::vector<Tape::Term>& terms, int slot, double weight);

    /*!
     * Add weight times the data point to the sum of data terms.
     */
    static void add_data_term(std::vector<Tape::DataTerm>& data_terms, int point, double weight);

    /*!
     * Return the set of boundaries, as a bit set of location indices, that the
     * face i of component axis lies beyond.
     */
    unsigned int getRegion(unsigned int axis, const hier::Index<NDIM>& i) const;

    /*!
     * Return the slot of the face i of component axis in the patch data, or
     * NOT_COMPUTABLE if the face is not in the patch data.
     */
    int getFaceSlot(unsigned int axis, const hier::Index<NDIM>& i);

    /*!
     * Return the slot of the value of the face i of component axis, and add
     * the assignments that compute it to the tape.
     */
    int getValue(unsigned int axis, const hier::Index<NDIM>& i);

    /*!
     * Return the slot of the value that the boundary loc gives to the face i,
     * and add the assignments that compute it to the tape.
     */
    int getValueAcross(unsigned int loc, unsigned int axis, const hier::Index<NDIM>& i);

    /*!
     * Return the slot of the value of the face i of a cell in region: the value
     * that the boundary loc gives to it if it lies beyond the same boundaries as
     * the cell, and the value it has otherwise.
     */
    int getFieldValue(unsigned int loc, unsigned int region, unsigned int axis, const hier::Index<NDIM>& i);

    /*!
     * Set value to the value that the boundary loc gives to the face i of the
     * component normal to it, which makes the divergence of the adjacent ghost
     * cell vanish, and return false if it is not computable.
     */
    bool setNormalValue(unsigned int loc, unsigned int axis, const hier::Index<NDIM>& i, Combination& value);

    /*!
     * Set value to the value that the boundary loc gives to the face i of a
     * component tangential to it, according to the boundary condition, and
     * return false if it is not computable.
     */
    bool setTangentialValue(unsigned int loc, unsigned int axis, const hier::Index<NDIM>& i, Combination& value);

    /*!
     * Add weight times the normal velocity at the face f of the boundary loc
     * to value, and return false if it is not computable.
     */
    bool addNormalVelocity(Combination& value,
                           const hier::Index<NDIM>& f,
                           unsigned int loc,
                           unsigned int axis,
                           unsigned int region,
                           double weight);

    /*!
     * Return how component axis is determined across the boundary loc at the
     * boundary face index.
     */
    TangentialRule getRule(unsigned int axis, unsigned int loc, const hier::Index<NDIM>& index) const;

    /*!
     * Return the index in the tape of the data point of component axis at the
     * boundary face index of the boundary loc, adding it if necessary.
     */
    int getDataPoint(unsigned int axis, unsigned int loc, const hier::Index<NDIM>& index, bool divide_by_mu);

    /*!
     * Evaluate the coefficients a and b of the boundary conditions.
     */
    void readRules();

    /*!
     * Add the assignment of value to the slot target to the tape.
     */
    void finish(Combination& value, int target);

    /*!
     * Group the data points of the tape into batches that one call to a
     * coefficient object evaluates.
     */
    void groupDataPoints();

    /*!
     * Return a new temporary slot.
     */
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

    std::map<std::pair<unsigned int, unsigned int>, Rules> d_rules;
    std::map<std::array<int, NDIM + 1>, int> d_face_ids;
    std::map<std::array<int, NDIM + 2>, int> d_point_ids;
    std::map<Key, int> d_memo;
    std::shared_ptr<Tape> d_tape;
};

/////////////////////////////// TAPE BUILDER /////////////////////////////////

INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::TapeBuilder(const INSStaggeredHierarchyIntegrator& fluid_solver,
                                                               traction_stencil::PhysicalDomainCache& domain_cache,
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
      d_tape(std::make_shared<Tape>())
{
    constexpr const char* method = "INSStaggeredDivergenceFreePhysBdryOp";
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
    const BoxArray<NDIM>& domain = traction_stencil::get_physical_domain(domain_cache, hierarchy, pgeom->getRatio());
    Box<NDIM> domain_box = domain[0];
    for (int n = 1; n < domain.size(); ++n)
    {
        domain_box += domain[n];
    }
    BoxList<NDIM> uncovered(domain_box);
    uncovered.removeIntersections(BoxList<NDIM>(domain));
    if (uncovered.getNumberOfItems() != 0)
    {
        TBOX_ERROR(method << ": the physical domain must be a single rectangular box.\n");
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
            TBOX_ERROR(method << ": the physical domain is " << d_upper[d] - d_lower[d] + 1 << " cells wide along axis "
                              << d << ", but the ghost cell width of the velocity is " << data.getGhostCellWidth().max()
                              << ". The domain must be at least as wide as the ghost cell width.\n");
        }
    }
}

std::shared_ptr<INSStaggeredDivergenceFreePhysBdryOp::Tape>
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::build()
{
    readRules();
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (Box<NDIM>::Iterator b(d_fill_box[axis]); b; b++)
        {
            const hier::Index<NDIM>& i = b();
            if (getRegion(axis, i) != 0 && getValue(axis, i) == NOT_COMPUTABLE)
            {
                d_tape->not_computable_faces.push_back(-getFaceSlot(axis, i) - 1);
            }
        }
    }
    groupDataPoints();
    return d_tape;
}

unsigned int
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::getRegion(const unsigned int axis, const hier::Index<NDIM>& i) const
{
    unsigned int region = 0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d_periodic[d])
        {
            continue;
        }
        if (i(d) < d_lower[d])
        {
            region |= 1u << location_index(d, false);
        }
        else if (i(d) > d_upper[d] + (d == axis ? 1 : 0))
        {
            region |= 1u << location_index(d, true);
        }
    }
    return region;
}

int
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::getFaceSlot(const unsigned int axis, const hier::Index<NDIM>& i)
{
    if (!d_data_box[axis].contains(i))
    {
        return NOT_COMPUTABLE;
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
    const int id = static_cast<int>(d_tape->faces.size());
    d_tape->faces.push_back({ axis, d_data_box[axis].offset(i) });
    d_face_ids[key] = id;
    return -id - 1;
}

int
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::newTemporary()
{
    return d_tape->num_temporaries++;
}

void
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::add_term(std::vector<Tape::Term>& terms,
                                                            const int slot,
                                                            const double weight)
{
    for (auto& term : terms)
    {
        if (term.slot == slot)
        {
            term.weight += weight;
            return;
        }
    }
    terms.push_back({ slot, weight });
}

void
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::add_data_term(std::vector<Tape::DataTerm>& data_terms,
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

int
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::getValue(const unsigned int axis, const hier::Index<NDIM>& i)
{
    if (!d_data_box[axis].contains(i))
    {
        return NOT_COMPUTABLE;
    }
    const unsigned int region = getRegion(axis, i);
    if (region == 0)
    {
        return getFaceSlot(axis, i);
    }
    Key key;
    key[0] = 0;
    key[1] = static_cast<int>(axis);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        key[d + 2] = i(d);
    }
    const auto it = d_memo.find(key);
    if (it != d_memo.end())
    {
        return it->second;
    }

    // Faces that are stored are assigned in the patch data, and faces that are
    // only needed to compute others are assigned temporary slots.
    const int target = d_fill_box[axis].contains(i) ? getFaceSlot(axis, i) : newTemporary();
    const int n = count_boundaries(region);
    Combination value;
    bool computable = true;
    for (unsigned int loc = 0; loc < 2 * NDIM && computable; ++loc)
    {
        if ((region & (1u << loc)) == 0)
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
            const int w = getValueAcross(loc, axis, i);
            computable = w != NOT_COMPUTABLE;
            if (computable)
            {
                add_term(value.terms, w, 1.0 / n);
            }
        }
    }
    int result = NOT_COMPUTABLE;
    if (computable)
    {
        finish(value, target);
        result = target;
    }
    d_memo[key] = result;
    return result;
}

int
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::getValueAcross(const unsigned int loc,
                                                                  const unsigned int axis,
                                                                  const hier::Index<NDIM>& i)
{
    const unsigned int region = getRegion(axis, i);
    if (count_boundaries(region) == 1)
    {
        return getValue(axis, i);
    }
    if (!d_data_box[axis].contains(i))
    {
        return NOT_COMPUTABLE;
    }
    Key key;
    key[0] = 1 + static_cast<int>(loc);
    key[1] = static_cast<int>(axis);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        key[d + 2] = i(d);
    }
    const auto it = d_memo.find(key);
    if (it != d_memo.end())
    {
        return it->second;
    }
    const int target = newTemporary();
    Combination value;
    const bool computable =
        (loc / 2 == axis) ? setNormalValue(loc, axis, i, value) : setTangentialValue(loc, axis, i, value);
    int result = NOT_COMPUTABLE;
    if (computable)
    {
        finish(value, target);
        result = target;
    }
    d_memo[key] = result;
    return result;
}

int
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::getFieldValue(const unsigned int loc,
                                                                 const unsigned int region,
                                                                 const unsigned int axis,
                                                                 const hier::Index<NDIM>& i)
{
    return getRegion(axis, i) == region ? getValueAcross(loc, axis, i) : getValue(axis, i);
}

bool
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::setNormalValue(const unsigned int loc,
                                                                  const unsigned int axis,
                                                                  const hier::Index<NDIM>& i,
                                                                  Combination& value)
{
    const bool is_upper = loc % 2 == 1;
    const unsigned int region = getRegion(axis, i);
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
    const int s_inner = getFieldValue(loc, region, axis, inner);
    if (s_inner == NOT_COMPUTABLE)
    {
        return false;
    }
    add_term(value.terms, s_inner, 1.0);
    for (unsigned int t = 0; t < NDIM; ++t)
    {
        if (t == axis)
        {
            continue;
        }
        hier::Index<NDIM> cell_hi = cell;
        cell_hi(t) += 1;
        const int s_hi = getFieldValue(loc, region, t, cell_hi);
        const int s_lo = getFieldValue(loc, region, t, cell);
        if (s_hi == NOT_COMPUTABLE || s_lo == NOT_COMPUTABLE)
        {
            return false;
        }
        const double weight = (is_upper ? -1.0 : 1.0) * d_dx[axis] / d_dx[t];
        add_term(value.terms, s_hi, weight);
        add_term(value.terms, s_lo, -weight);
    }
    return true;
}

bool
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::setTangentialValue(const unsigned int loc,
                                                                      const unsigned int axis,
                                                                      const hier::Index<NDIM>& i,
                                                                      Combination& value)
{
    const unsigned int normal_axis = loc / 2;
    const bool is_upper = loc % 2 == 1;
    const unsigned int region = getRegion(axis, i);
    const int k = is_upper ? i(normal_axis) - d_upper[normal_axis] - 1 : d_lower[normal_axis] - 1 - i(normal_axis);
    hier::Index<NDIM> mirror = i;
    mirror(normal_axis) = is_upper ? d_upper[normal_axis] - k : d_lower[normal_axis] + k;
    const int s_mirror = getValue(axis, mirror);
    if (s_mirror == NOT_COMPUTABLE)
    {
        return false;
    }
    hier::Index<NDIM> foot = i;
    foot(normal_axis) = is_upper ? d_upper[normal_axis] + 1 : d_lower[normal_axis];
    const double sgn = is_upper ? +1.0 : -1.0;
    const double gamma_weight = (2.0 * k + 1.0) * d_dx[normal_axis] * sgn;
    switch (getRule(axis, loc, foot))
    {
    case TangentialRule::PRESCRIBED:
        add_term(value.terms, s_mirror, -1.0);
        add_data_term(value.data_terms, getDataPoint(axis, loc, foot, false), 2.0);
        break;
    case TangentialRule::PSEUDO_TRACTION:
        add_term(value.terms, s_mirror, 1.0);
        add_data_term(value.data_terms, getDataPoint(axis, loc, foot, true), gamma_weight);
        break;
    case TangentialRule::TRACTION:
    {
        add_term(value.terms, s_mirror, 1.0);
        add_data_term(value.data_terms, getDataPoint(axis, loc, foot, true), gamma_weight);
        hier::Index<NDIM> foot_lower = foot;
        foot_lower(axis) -= 1;
        const double weight = gamma_weight / d_dx[axis];
        if (!addNormalVelocity(value, foot, loc, axis, region, -weight) ||
            !addNormalVelocity(value, foot_lower, loc, axis, region, +weight))
        {
            return false;
        }
        break;
    }
    }
    return true;
}

bool
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::addNormalVelocity(Combination& value,
                                                                     const hier::Index<NDIM>& f,
                                                                     const unsigned int loc,
                                                                     const unsigned int axis,
                                                                     const unsigned int region,
                                                                     const double weight)
{
    const unsigned int normal_axis = loc / 2;
    const unsigned int other_region = region & ~(1u << loc);
    const unsigned int f_region = getRegion(normal_axis, f);
    const unsigned int extra = f_region & ~other_region;

    // A face that is not beyond a corner has its value. A face one cell beyond a
    // corner, across a boundary that is not in the region, has the value that
    // the corner rule gives: twice the velocity that the adjacent boundary
    // prescribes minus the value at the nearest face on B, or the linear
    // extrapolation from the faces on B if the adjacent boundary does not
    // prescribe the normal velocity.
    if (extra == 0)
    {
        const int s = getValue(normal_axis, f);
        if (s == NOT_COMPUTABLE)
        {
            return false;
        }
        add_term(value.terms, s, weight);
        return true;
    }
    const unsigned int adj_loc_lower = location_index(axis, false);
    const unsigned int adj_loc_upper = location_index(axis, true);
    if (extra != (1u << adj_loc_lower) && extra != (1u << adj_loc_upper))
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp: the normal velocity face "
                   << f << " is not one cell beyond a corner of the boundary " << loc << ".\n");
    }
    const int step = extra == (1u << adj_loc_lower) ? +1 : -1;

    // Linear extrapolation uses the faces on B that are two cells inward from f along axis; a boundary segment that is
    // one cell long has no second face.
    const int j_next = f(axis) + 2 * step;
    const bool can_extrapolate = j_next >= d_lower[axis] && j_next <= d_upper[axis];
    traction_stencil::NormalVelocityStencil stencil =
        traction_stencil::get_corner_normal_velocity_stencil(f, normal_axis, axis, step, can_extrapolate);
    hier::Index<NDIM> adj_foot = f;
    adj_foot(axis) = stencil.adjacent_location_index % 2 == 1 ? d_upper[axis] + 1 : d_lower[axis];
    if (getRule(normal_axis, stencil.adjacent_location_index, adj_foot) == TangentialRule::PRESCRIBED)
    {
        stencil.useAdjacentBoundaryValue();
    }
    for (int n = 0; n < 2; ++n)
    {
        if (stencil.weight[n] == 0.0)
        {
            continue;
        }
        const int s = getValue(normal_axis, stencil.idx[n]);
        if (s == NOT_COMPUTABLE)
        {
            return false;
        }
        add_term(value.terms, s, weight * stencil.weight[n]);
    }
    if (stencil.uses_adjacent_boundary_value)
    {
        add_data_term(value.data_terms,
                      getDataPoint(normal_axis, stencil.adjacent_location_index, adj_foot, false),
                      2.0 * weight);
    }
    return true;
}

void
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::readRules()
{
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs = d_fluid_solver.getPhysicalBoundaryConditions();
    Box<NDIM> cell_box = d_data_cell_box;
    cell_box.grow(1);

    // The coefficients of a component are evaluated at its faces, with the geometry of the patch shifted by half a
    // cell.
    std::array<std::unique_ptr<traction_stencil::ShiftedPatchGeometry>, NDIM> shifted_geometry;
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
                    TBOX_ERROR(
                        "INSStaggeredDivergenceFreePhysBdryOp: the fluid solver has no velocity boundary "
                        "conditions for component "
                        << axis << ".\n");
                }
                Rules rules;
                rules.box = SideGeometry<NDIM>::toSideBox(cell_box, axis);
                const int plane = is_upper ? d_upper[normal_axis] + 1 : d_lower[normal_axis];
                rules.box.lower()(normal_axis) = plane;
                rules.box.upper()(normal_axis) = plane;
                rules.acoef_data = new ArrayData<NDIM, double>(rules.box, 1);
                rules.bcoef_data = new ArrayData<NDIM, double>(rules.box, 1);
                Pointer<ArrayData<NDIM, double>> gcoef_data;
                if (!shifted_geometry[axis])
                {
                    std::array<double, NDIM> shift;
                    shift.fill(0.0);
                    shift[axis] = -0.5 * d_dx[axis];
                    shifted_geometry[axis] = std::make_unique<traction_stencil::ShiftedPatchGeometry>(d_patch, shift);
                }
                Box<NDIM> bdry_cell_box = rules.box;
                if (!is_upper)
                {
                    bdry_cell_box.shift(normal_axis, -1);
                }
                const BoundaryBox<NDIM> bdry_box(bdry_cell_box, 1, loc);
                {
                    const traction_stencil::ShiftedPatchGeometry::Scope scope(d_patch, *shifted_geometry[axis]);
                    bc_coefs[axis]->setBcCoefs(
                        rules.acoef_data, rules.bcoef_data, gcoef_data, d_var, d_patch, bdry_box, d_fill_time);
                }
                d_rules[std::make_pair(axis, loc)] = rules;
            }
        }
    }
}

TangentialRule
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::getRule(const unsigned int axis,
                                                           const unsigned int loc,
                                                           const hier::Index<NDIM>& index) const
{
    const auto it = d_rules.find(std::make_pair(axis, loc));
    if (it == d_rules.end() || !it->second.box.contains(index))
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp: the boundary condition of component "
                   << axis << " on the boundary " << loc << " is not known at the boundary face " << index << ".\n");
    }
    const double alpha = (*it->second.acoef_data)(index, 0);
    const double beta = (*it->second.bcoef_data)(index, 0);
    const bool velocity_bc = (alpha == 1.0 && beta == 0.0);
    const bool traction_bc = (alpha == 0.0 && beta == 1.0);
    if (!velocity_bc && !traction_bc)
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp: the boundary condition of component "
                   << axis << " on the boundary " << loc << " at the boundary face " << index
                   << " must have (a, b) = (1, 0) or (a, b) = (0, 1), but (a, b) = (" << alpha << ", " << beta
                   << ").\n");
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
        TBOX_ERROR(
            "INSStaggeredDivergenceFreePhysBdryOp: unrecognized or unsupported traction boundary condition "
            "type: "
            << enum_to_string<TractionBcType>(d_traction_bc_type) << "\n");
    }
    return TangentialRule::TRACTION;
}

int
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::getDataPoint(const unsigned int axis,
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
    const int id = static_cast<int>(d_tape->points.size());
    d_tape->points.push_back({ axis, loc, index, divide_by_mu });
    d_point_ids[key] = id;
    return id;
}

void
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::finish(Combination& value, const int target)
{
    Tape::Entry entry;
    entry.target = target;
    entry.first_term = d_tape->terms.size();
    entry.num_terms = value.terms.size();
    entry.first_data_term = d_tape->data_terms.size();
    entry.num_data_terms = value.data_terms.size();
    d_tape->terms.insert(d_tape->terms.end(), value.terms.begin(), value.terms.end());
    d_tape->data_terms.insert(d_tape->data_terms.end(), value.data_terms.begin(), value.data_terms.end());
    d_tape->entries.push_back(entry);
}

void
INSStaggeredDivergenceFreePhysBdryOp::TapeBuilder::groupDataPoints()
{
    std::map<std::pair<unsigned int, unsigned int>, std::size_t> batch_ids;
    for (int p = 0; p < static_cast<int>(d_tape->points.size()); ++p)
    {
        const Tape::DataPoint& point = d_tape->points[p];
        const auto key = std::make_pair(point.axis, point.location_index);
        auto it = batch_ids.find(key);
        if (it == batch_ids.end())
        {
            it = batch_ids.insert(std::make_pair(key, d_tape->batches.size())).first;
            d_tape->batches.push_back({ point.axis, point.location_index, Box<NDIM>(point.index, point.index), {} });
        }
        Tape::Batch& batch = d_tape->batches[it->second];
        batch.box += Box<NDIM>(point.index, point.index);
        batch.points.push_back(p);
    }
}

/////////////////////////////// TAPE /////////////////////////////////////////

std::vector<double>
INSStaggeredDivergenceFreePhysBdryOp::Tape::evaluateDataPoints(Patch<NDIM>& patch,
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
    std::array<std::unique_ptr<traction_stencil::ShiftedPatchGeometry>, NDIM> shifted_geometry;
    for (const auto& batch : batches)
    {
        if (!shifted_geometry[batch.axis])
        {
            std::array<double, NDIM> shift;
            shift.fill(0.0);
            shift[batch.axis] = -0.5 * pgeom->getDx()[batch.axis];
            shifted_geometry[batch.axis] = std::make_unique<traction_stencil::ShiftedPatchGeometry>(patch, shift);
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
            const traction_stencil::ShiftedPatchGeometry::Scope scope(patch, *shifted_geometry[batch.axis]);
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
INSStaggeredDivergenceFreePhysBdryOp::Tape::apply(const std::array<double*, NDIM>& u,
                                                  const std::vector<double>& point_values) const
{
    std::vector<double> temporaries(num_temporaries, 0.0);
    const auto value = [&](const int slot) -> double&
    {
        if (slot >= 0)
        {
            return temporaries[slot];
        }
        const Face& face = faces[-slot - 1];
        return u[face.axis][face.offset];
    };
    for (const auto& entry : entries)
    {
        double sum = 0.0;
        for (std::size_t n = entry.first_term; n < entry.first_term + entry.num_terms; ++n)
        {
            sum += terms[n].weight * value(terms[n].slot);
        }
        for (std::size_t n = entry.first_data_term; n < entry.first_data_term + entry.num_data_terms; ++n)
        {
            sum += data_terms[n].weight * point_values[data_terms[n].point];
        }
        value(entry.target) = sum;
    }
    for (const int id : not_computable_faces)
    {
        u[faces[id].axis][faces[id].offset] = std::numeric_limits<double>::quiet_NaN();
    }
}

void
INSStaggeredDivergenceFreePhysBdryOp::Tape::applyTranspose(const std::array<double*, NDIM>& u) const
{
    std::vector<double> temporaries(num_temporaries, 0.0);
    const auto value = [&](const int slot) -> double&
    {
        if (slot >= 0)
        {
            return temporaries[slot];
        }
        const Face& face = faces[-slot - 1];
        return u[face.axis][face.offset];
    };
    for (auto it = entries.rbegin(); it != entries.rend(); ++it)
    {
        const double accumulated = value(it->target);
        value(it->target) = 0.0;
        for (std::size_t n = it->first_term; n < it->first_term + it->num_terms; ++n)
        {
            value(terms[n].slot) += terms[n].weight * accumulated;
        }
    }
}

int
INSStaggeredDivergenceFreePhysBdryOp::Tape::findNonzeroNotComputableFace(const std::array<double*, NDIM>& u) const
{
    for (const int id : not_computable_faces)
    {
        if (u[faces[id].axis][faces[id].offset] != 0.0)
        {
            return id;
        }
    }
    return -1;
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
    d_tape_cache.clear();
    return;
} // clearCache

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
        const std::shared_ptr<const Tape> tape = getTape(patch, patch_data_idx, ghost_width_to_fill, fill_time);
        Pointer<SideData<NDIM, double>> data = patch.getPatchData(patch_data_idx);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<Variable<NDIM>> var;
        var_db->mapIndexToVariable(patch_data_idx, var);
        const std::vector<double> point_values = tape->evaluateDataPoints(
            patch, var, d_fluid_solver->getPhysicalBoundaryConditions(), fill_time, mu, d_homogeneous_bc);
        std::array<double*, NDIM> u;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            u[axis] = data->getPointer(axis, 0);
        }
        tape->apply(u, point_values);
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
        const std::shared_ptr<const Tape> tape = getTape(patch, patch_data_idx, ghost_width_to_fill, fill_time);
        Pointer<SideData<NDIM, double>> data = patch.getPatchData(patch_data_idx);
        std::array<double*, NDIM> u;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            u[axis] = data->getPointer(axis, 0);
        }
        const int id = tape->findNonzeroNotComputableFace(u);
        if (id >= 0)
        {
            const Tape::Face& face = tape->faces[id];
            TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  patch " << patch.getPatchNumber() << " (box " << patch.getBox() << ") has the value "
                       << u[face.axis][face.offset] << " at the face "
                       << data->getArrayData(face.axis).getBox().index(face.offset) << " with axis " << face.axis
                       << ", where the extension needs values outside the patch data and cannot be transposed.\n");
        }
        tape->applyTranspose(u);
    }
    return;
} // accumulateFromPhysicalBoundaryData

/////////////////////////////// PRIVATE //////////////////////////////////////

std::shared_ptr<const INSStaggeredDivergenceFreePhysBdryOp::Tape>
INSStaggeredDivergenceFreePhysBdryOp::getTape(Patch<NDIM>& patch,
                                              const int patch_data_idx,
                                              const IntVector<NDIM>& ghost_width_to_fill,
                                              const double fill_time)
{
    Pointer<SideData<NDIM, double>> data = patch.getPatchData(patch_data_idx);
    if (!data)
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp:\n"
                   << "  patch data index " << patch_data_idx
                   << " does not correspond to a side-centered double precision variable.\n");
    }
    if (data->getDepth() != 1)
    {
        TBOX_ERROR("INSStaggeredDivergenceFreePhysBdryOp:\n"
                   << "  the depth of patch data index " << patch_data_idx << " is " << data->getDepth()
                   << " but a velocity has depth 1.\n");
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
        const auto it = d_tape_cache.find(key);
        if (it != d_tape_cache.end())
        {
            return it->second;
        }
    }
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM>> var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    TapeBuilder builder(*d_fluid_solver, d_physical_domain, patch, var, *data, gcw_to_fill, fill_time);
    std::shared_ptr<const Tape> tape = builder.build();
    if (in_hierarchy)
    {
        d_tape_cache[key] = tape;
    }
    return tape;
} // getTape

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
