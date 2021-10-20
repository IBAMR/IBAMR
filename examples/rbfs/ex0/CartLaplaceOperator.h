// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_CartLaplaceOperator
#define included_IBTK_CartLaplaceOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/PETScLinearAugmentedOperator.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "VariableFillPattern.h"
#include "tbox/Pointer.h"

#include "libmesh/boundary_mesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/node.h"

#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
struct IndexList
{
public:
    IndexList(const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch, const SAMRAI::pdat::CellIndex<NDIM>& idx)
        : d_idx(idx), d_patch(patch)
    {
        const SAMRAI::hier::Box<NDIM>& box = patch->getBox();
        const SAMRAI::hier::Index<NDIM>& idx_low = box.lower();
        const SAMRAI::hier::Index<NDIM>& idx_up = box.upper();
        int num_x = idx_up(0) - idx_low(0) + 1;
        d_global_idx = idx(0) - idx_low(0) + num_x * (idx(1) - idx_low(1) + 1);
#if (NDIM == 3)
        int num_y = idx_up(1) - idx_low(1) + 1;
        d_global_idx += num_x * num_y * (idx(2) - idx_low(2));
#endif
    }

    bool operator<(const IndexList& b) const
    {
        bool less_than_b = false;
        if (d_patch->getPatchNumber() < b.d_patch->getPatchNumber())
            less_than_b = true;
        else if (b.d_patch->getPatchNumber() == d_patch->getPatchNumber() && d_global_idx < b.d_global_idx)
        {
            less_than_b = true;
        }
        return less_than_b;
    }
    int d_global_idx = -1;
    SAMRAI::pdat::CellIndex<NDIM> d_idx;
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > d_patch;
};

struct UPoint
{
public:
    UPoint(const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch, const SAMRAI::pdat::CellIndex<NDIM>& idx)
        : d_idx(idx)
    {
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        const double* const xlow = pgeom->getXLower();
        const SAMRAI::hier::Index<NDIM>& idx_low = patch->getBox().lower();
        for (unsigned int d = 0; d < NDIM; ++d)
            d_pt[d] = xlow[d] + dx[d] * (static_cast<double>(d_idx(d) - idx_low(d)) + 0.5);
    };

    UPoint(const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch, libMesh::Node* node) : d_node(node)
    {
        for (unsigned int d = 0; d < NDIM; ++d) d_pt[d] = (*d_node)(d);
    };

    UPoint() : d_empty(true)
    {
        // intentionally blank
    }

    UPoint(const std::vector<double>& pt) : d_empty(true)
    {
        for (unsigned int d = 0; d < NDIM; ++d) d_pt[d] = pt[d];
    }

    double dist(const IBTK::VectorNd& x) const
    {
        return (d_pt - x).norm();
    }

    double dist(const UPoint& pt) const
    {
        return (d_pt - pt.getVec()).norm();
    }

    double operator()(const size_t i) const
    {
        return d_pt(i);
    }

    double operator[](const size_t i) const
    {
        return d_pt[i];
    }

    friend std::ostream& operator<<(std::ostream& out, const UPoint& pt)
    {
        out << "   location: " << pt.d_pt.transpose() << "\n";
        if (pt.isNode())
            out << "   node id: " << pt.d_node->id() << "\n";
        else if (!pt.isEmpty())
            out << "   idx:     " << pt.d_idx << "\n";
        else
            out << "   pt is neither node nor index\n";
        return out;
    }

    bool isEmpty() const
    {
        return d_empty;
    }

    bool isNode() const
    {
        return d_node != nullptr;
    }

    const libMesh::Node* const getNode() const
    {
        if (!isNode()) TBOX_ERROR("Not at a node\n");
        return d_node;
    }

    const SAMRAI::pdat::CellIndex<NDIM>& getIndex() const
    {
        if (isNode()) TBOX_ERROR("At at node\n");
        if (d_empty) TBOX_ERROR("Not a point\n");
        return d_idx;
    }

    const VectorNd& getVec() const
    {
        return d_pt;
    }

private:
    IBTK::VectorNd d_pt;
    libMesh::Node* d_node = nullptr;
    SAMRAI::pdat::CellIndex<NDIM> d_idx;
    bool d_empty = false;
};

/*!
 * \brief Class CartLaplaceOperator is a concrete LaplaceOperator which implements
 * a globally second-order accurate cell-centered finite difference
 * discretization of a scalar elliptic operator of the form \f$ L = C I + \nabla
 * \cdot D \nabla\f$.
 */
class CartLaplaceOperator : public PETScLinearAugmentedOperator
{
public:
    /*!
     * \brief Constructor for class LaplaceOperator initializes the operator
     * coefficients and boundary conditions to default values.
     */
    CartLaplaceOperator(const std::string& object_name,
                        libMesh::BoundaryMesh* bdry_mesh,
                        const libMesh::DofMap* dof_map,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~CartLaplaceOperator();

    /*!
     * \brief Set the level set
     */
    inline void setLS(int ls_idx)
    {
        d_ls_idx = ls_idx;
    }

    void setupBeforeApply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                          SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y);

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Compute y=Ax.
     *
     * Before calling this function, the form of the vectors x and y should be
     * set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * Thus, the user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must
     * have enough ghost cells for computation of Ax).
     *
     * \note In general, the vectors x and y \em cannot be the same.
     *
     * Upon return from this function, the y vector will contain the result of
     * the application of A to x.
     *
     * initializeOperatorState must be called prior to any calls to
     * applyOperator.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=Ax
     */
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Compute hierarchy-dependent data required for computing y=Ax (and
     * y=A'x).
     *
     * \param in input vector
     * \param out output vector
     *
     * \see KrylovLinearSolver::initializeSolverState
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    /*!
     * \brief Remove all hierarchy-dependent data computed by
     * initializeOperatorState().
     *
     * Remove all hierarchy-dependent data set by initializeOperatorState().  It
     * is safe to call deallocateOperatorState() even if the state is already
     * deallocated.
     *
     * \see initializeOperatorState
     * \see KrylovLinearSolver::deallocateSolverState
     */
    void deallocateOperatorState() override;

    //\}

    /*!
     * Debugging functions
     */
    void sortLagDOFsToCells();
    void printPtMap(std::ostream& os);
    inline void setPatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy)
    {
        d_hierarchy = patch_hierarchy;
    }

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CartLaplaceOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CartLaplaceOperator(const LaplaceOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LaplaceOperator& operator=(const LaplaceOperator& that) = delete;

    void applyToLagDOFs(int x_idx, int y_idx);

    double formRBFDer();

    double getSolVal(const UPoint& pt, const SAMRAI::pdat::CellData<NDIM, double>& Q_data) const;
    void setSolVal(double q, const UPoint& pt, SAMRAI::pdat::CellData<NDIM, double>& Q_data) const;

    static unsigned int s_num_ghost_cells;
    double d_eps = std::numeric_limits<double>::quiet_NaN();

    // Operator parameters.
    int d_ncomp = 0;

    // Cached communications operators.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_fill_pattern;
    std::vector<HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_transaction_comps;
    SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> d_hier_bdry_fill, d_no_fill;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_x, d_b;

    int d_ls_idx = IBTK::invalid_index;
    double d_dist_to_bdry = std::numeric_limits<double>::quiet_NaN();

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Lag structure info
    libMesh::BoundaryMesh* d_mesh;
    std::vector<std::vector<libMesh::Node*> > d_idx_node_vec, d_idx_node_ghost_vec;
    std::vector<std::vector<UPoint> > d_base_pt_vec;
    std::vector<std::vector<std::vector<UPoint> > > d_pair_pt_vec;

    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    const libMesh::DofMap* d_dof_map = nullptr;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CartLaplaceOperator
