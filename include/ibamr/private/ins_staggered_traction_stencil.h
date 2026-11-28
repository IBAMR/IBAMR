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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_ins_staggered_traction_stencil
#define included_IBAMR_ins_staggered_traction_stencil

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/StokesBcCoefStrategy.h>

#include <ibtk/ibtk_enums.h>

#include <tbox/Pointer.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideIndex.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace detail
{
class ScopedPatchCopy;

/*!
 * \brief Cache of the periodic physical domain and the patch copies of the TRACTION correction, keyed by value.
 */
class TractionPatchCache
{
public:
    /*!
     * \brief The largest number of patch copies, and of pairs of coefficient arrays, held.
     */
    static constexpr std::size_t max_cached_entries = 256;

    /*!
     * \brief Return the physical domain of the grid geometry of \a hierarchy, refined to the level of \a patch,
     * together with its periodic images. It is built once for each refinement ratio, periodic shift, and coarse
     * physical domain. The reference remains valid for the life of the cache.
     */
    const SAMRAI::hier::BoxArray<NDIM>&
    getPhysicalDomain(const SAMRAI::hier::Patch<NDIM>& patch,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy);

    /*!
     * \brief Return the coefficient arrays of depth one on \a box for evaluating the dependence of the boundary data on
     * the velocity. Each call sets every value of the arrays to the largest double, the value of a new array of an
     * IBSAMRAI2 debug build. The arrays are allocated once for each box and reused, as the patch copies are.
     */
    std::pair<SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>>
    getCoefficientData(const SAMRAI::hier::Box<NDIM>& box);

private:
    friend class ScopedPatchCopy;

    /*!
     * A map from the indices of a box to values, with at most max_cached_entries entries. Each find() or insert()
     * counts as a use of its entry. Inserting into a full map erases the less recently used half of the entries. The
     * use counts are values.
     */
    template <class Value>
    class UseMap
    {
    public:
        using Key = std::array<int, 2 * NDIM>;

        /*!
         * Return the value of \a key, or nullptr if there is none.
         */
        Value* find(const Key& key);

        /*!
         * Add \a value for \a key, which must have no entry, and return it.
         */
        Value& insert(const Key& key, Value value);

        /*!
         * Erase the entry of \a key.
         */
        void erase(const Key& key);

    private:
        struct Entry
        {
            Value value;
            std::uint64_t last_use;
        };

        std::map<Key, Entry> d_entries;
        std::uint64_t d_use_count = 0;
    };

    /*!
     * Return a copy of \a patch that shares its patch data and its geometry. The copy is allocated once for each box,
     * and each call sets its geometry and patch data to those of \a patch. A copy is reused only for an equal box and
     * patch descriptor.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> getPatchCopy(const SAMRAI::hier::Patch<NDIM>& patch);

    /*!
     * A cached domain and the values that determine it.
     */
    struct Domain
    {
        bool matches(const SAMRAI::hier::IntVector<NDIM>& ratio,
                     const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                     const SAMRAI::hier::BoxArray<NDIM>& coarse_domain) const;

        SAMRAI::hier::IntVector<NDIM> ratio;
        SAMRAI::hier::IntVector<NDIM> periodic_shift;
        SAMRAI::hier::BoxArray<NDIM> coarse_domain;
        SAMRAI::hier::BoxArray<NDIM> domain;
    };

    // References to the elements of a deque stay valid as it grows.
    std::deque<Domain> d_domains;

    // The key is the lower and upper indices of the box. IBSAMRAI2's IntVector operator< is componentwise, not a strict
    // weak order, so a std::array of the indices serves as the key.
    UseMap<SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>> d_copies;
    UseMap<std::pair<SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>>>
        d_coefficient_data;
};

/*!
 * \brief A copy of a patch from a TractionPatchCache for the duration of the object. The copy shares the patch data
 * of the patch, and replacing the geometry of the copy leaves the patch unchanged. The destructor releases the patch
 * data and the geometry of the copy, so that their lifetime does not depend on the copy. The constructor does not set
 * the patch number, level number, or in-hierarchy flag of the copy.
 */
class ScopedPatchCopy
{
public:
    /*!
     * \brief Constructor.
     */
    ScopedPatchCopy(TractionPatchCache& cache, const SAMRAI::hier::Patch<NDIM>& patch);

    /*!
     * \brief Destructor.
     */
    ~ScopedPatchCopy();

    ScopedPatchCopy(const ScopedPatchCopy&) = delete;
    ScopedPatchCopy& operator=(const ScopedPatchCopy&) = delete;

    /*!
     * \brief Return the copy.
     */
    SAMRAI::hier::Patch<NDIM>& get() const;

private:
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> d_copy;
};

/*!
 * \brief A copy of a patch for evaluating the boundary condition coefficients of the normal velocity component
 * \a bdry_normal_axis at the faces of the tangential component \a tangential_axis. The original patch is unchanged.
 */
class NormalVelocityPatch
{
public:
    /*!
     * \brief Constructor. The copy shares the patch data of \a patch, and the copy comes from \a cache.
     */
    NormalVelocityPatch(TractionPatchCache& cache,
                        const SAMRAI::hier::Patch<NDIM>& patch,
                        unsigned int bdry_normal_axis,
                        unsigned int tangential_axis);

    /*!
     * \brief Return the box of the patch.
     */
    const SAMRAI::hier::Box<NDIM>& getBox() const;

    /*!
     * \brief Return the axis normal to the physical boundary.
     */
    unsigned int getBdryNormalAxis() const;

    /*!
     * \brief Return the axis of the tangential velocity component.
     */
    unsigned int getTangentialAxis() const;

    /*!
     * \brief Return the grid spacing along the tangential axis.
     */
    double getTangentialDx() const;

    /*!
     * \brief Set the coefficients of \a normal_bc_coef on the boundary box \a bdry_box at the faces of the tangential
     * component, using the geometry of the copy shifted by half a cell along the tangential axis and minus half a cell
     * along the axis normal to the boundary.
     */
    void setBcCoefs(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>& acoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>& bcoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>& gcoef_data,
                    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                    double fill_time);

private:
    /*!
     * The copy of the patch and the axes of its normal and tangential velocity components.
     */
    ScopedPatchCopy d_patch;
    unsigned int d_bdry_normal_axis;
    unsigned int d_tangential_axis;
};

/*!
 * \brief Stencil weight[0]*u(idx[0]) + weight[1]*u(idx[1]) for the normal velocity in the TRACTION condition at a face
 * on a boundary. Beyond a corner, the weights are those of linear extrapolation, or of constant extrapolation from a
 * boundary segment of one face. If the normal velocity u_b is prescribed on the adjacent boundary
 * (adjacent_location_index), the caller calls useAdjacentBoundaryValue() and adds prescribed_weight*u_b.
 */
struct NormalVelocityStencil
{
    /*!
     * \brief The faces of the stencil.
     */
    std::array<SAMRAI::pdat::SideIndex<NDIM>, 2> idx;

    /*!
     * \brief The weights of the values at the faces.
     */
    std::array<double, 2> weight;

    /*!
     * \brief Whether the face lies beyond a corner of the physical boundary.
     */
    bool beyond_corner = false;

    /*!
     * \brief The location index of the adjacent boundary, if beyond_corner is true.
     */
    unsigned int adjacent_location_index = 0;

    /*!
     * \brief Whether useAdjacentBoundaryValue() has been called.
     */
    bool uses_adjacent_boundary_value = false;

    /*!
     * \brief The weight of the prescribed normal velocity u_b on the adjacent boundary, which
     * useAdjacentBoundaryValue() sets.
     */
    double prescribed_weight = 2.0;

    /*!
     * \brief Replace the weighted sum by -q_0, or by -2*q_0 + q_1/3 if \a bdry_interp_type is QUADRATIC and the
     * boundary segment has two faces, in which q_0 and q_1 are the values at the two nearest faces on the boundary, and
     * set prescribed_weight to 2 or 8/3. The member beyond_corner must be true.
     */
    void useAdjacentBoundaryValue(IBTK::BdryInterpType bdry_interp_type);
};

/*!
 * \brief Return the stencil of the normal velocity at the face \a i_face beyond a corner of the boundary normal to
 * \a bdry_normal_axis, where the boundary segment ends along \a tangential_axis. The nearest face on the boundary is
 * \a i_face + \a step along \a tangential_axis. The stencil also uses the second nearest face if \a can_extrapolate
 * is true.
 */
NormalVelocityStencil get_corner_normal_velocity_stencil(SAMRAI::hier::Index<NDIM> i_face,
                                                         unsigned int bdry_normal_axis,
                                                         unsigned int tangential_axis,
                                                         int step,
                                                         bool can_extrapolate);

/*!
 * \brief Return whether any face of \a acoef_data and \a bcoef_data has a TRACTION condition.
 */
bool has_traction_face(const SAMRAI::pdat::ArrayData<NDIM, double>& acoef_data,
                       const SAMRAI::pdat::ArrayData<NDIM, double>& bcoef_data);

/*!
 * \brief Correct the coefficients of the TRACTION condition mu*(du_tan/dx_norm + du_norm/dx_tan) = g on the tangential
 * velocity component \a tangential_axis at the boundary normal to \a bdry_normal_axis. The corrected coefficients are
 * (0, 1) and a value that includes the tangential derivative du_norm/dx_tan of the normal velocity. The coefficients
 * of faces with a velocity condition are unchanged.
 *
 * \a get_mu returns the viscosity at a face, given the traction datum g there. The function evaluates it only at faces
 * with a TRACTION condition.
 * \a patch_hierarchy must be set. The periodic physical domain and the copy of \a patch come from \a cache.
 *
 * Beyond a corner where the normal velocity is prescribed on the adjacent boundary, \a bdry_interp_type determines the
 * extrapolation of the normal velocity: see NormalVelocityStencil::useAdjacentBoundaryValue().
 */
void correct_traction_coefficients(SAMRAI::pdat::ArrayData<NDIM, double>& acoef_data,
                                   SAMRAI::pdat::ArrayData<NDIM, double>& bcoef_data,
                                   SAMRAI::pdat::ArrayData<NDIM, double>& gcoef_data,
                                   const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                   const SAMRAI::hier::Patch<NDIM>& patch,
                                   unsigned int bdry_normal_axis,
                                   bool bdry_is_lower,
                                   unsigned int tangential_axis,
                                   SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                   bool homogeneous_bc,
                                   const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& patch_hierarchy,
                                   TractionPatchCache& cache,
                                   const std::function<double(const SAMRAI::hier::Index<NDIM>&, double)>& get_mu,
                                   IBTK::BdryInterpType bdry_interp_type,
                                   double fill_time);

/*!
 * \brief Return the normal velocity on the face, farther from the boundary, of the ghost cell \a i_g outside the
 * physical boundary normal to \a bdry_normal_axis (the lower boundary if \a bdry_is_lower), such that the discrete
 * divergence of \a u_data vanishes in \a i_g. The function ignores \a u_data on that face.
 */
double get_divergence_free_ghost_value(const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                       const SAMRAI::hier::Index<NDIM>& i_g,
                                       unsigned int bdry_normal_axis,
                                       bool bdry_is_lower,
                                       const double* dx);

/*!
 * \brief Return the faces of the velocity on which get_divergence_free_ghost_value() depends at the ghost cell \a i_g,
 * each with its coefficient.
 */
StokesBcCoefStrategy::SideStencil get_divergence_free_ghost_value_stencil(const SAMRAI::hier::Index<NDIM>& i_g,
                                                                          unsigned int bdry_normal_axis,
                                                                          bool bdry_is_lower,
                                                                          const double* dx);

/*!
 * \brief Return whether any of \a u_bc_coefs is a StokesBcCoefStrategy with TRACTION conditions or, if
 * \a bdry_interp_type is QUADRATIC, PSEUDO_TRACTION conditions: the conditions for which
 * StaggeredStokesPhysicalBoundaryHelper::addNormalTractionViscousTerm() adds a term.
 */
bool has_normal_stress_term(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                            IBTK::BdryInterpType bdry_interp_type);

/*!
 * \brief Set the Robin coefficients of the velocity component \a axis on \a bdry_box of \a patch at the locations of
 * that component: the patch geometry is shifted by half a cell along \a axis for the call of \a bc_coef.
 */
void set_component_bc_coefs(SAMRAI::solv::RobinBcCoefStrategy<NDIM>& bc_coef,
                            unsigned int axis,
                            SAMRAI::hier::Patch<NDIM>& patch,
                            const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>& acoef_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>& bcoef_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>& gcoef_data,
                            double fill_time);

/*!
 * Return whether the TRACTION conditions of the velocity component with the Robin coefficient object \a comp_bc_coef
 * make its boundary data on \a bdry_box of \a patch at \a fill_time depend on the velocity: the component is tangential
 * to the box and a face that the fill sets has TRACTION coefficients. The copy of \a patch comes from \a cache.
 */
bool traction_data_depend_on_velocity(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* comp_bc_coef,
                                      unsigned int comp_idx,
                                      const SAMRAI::hier::Patch<NDIM>& patch,
                                      TractionPatchCache& cache,
                                      const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                                      double fill_time);

/*!
 * \brief Return the side indices of the normal velocities on which the inhomogeneous Robin coefficient of the TRACTION
 * condition depends at the location \a i, each with the derivative of the coefficient with respect to it. An index may
 * appear more than once: its derivative is the sum of its entries. The faces of the stencil lie in \a ghost_box, the
 * ghost box of the velocity data.
 *
 * Beyond a corner where the normal velocity is prescribed on the adjacent boundary, \a bdry_interp_type determines the
 * extrapolation of the normal velocity: see NormalVelocityStencil::useAdjacentBoundaryValue().
 */
StokesBcCoefStrategy::SideStencil get_traction_stencil(const SAMRAI::hier::Index<NDIM>& i,
                                                       bool bdry_is_lower,
                                                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                                       NormalVelocityPatch& normal_velocity_patch,
                                                       const SAMRAI::hier::BoxArray<NDIM>& domain,
                                                       const SAMRAI::hier::Box<NDIM>& ghost_box,
                                                       IBTK::BdryInterpType bdry_interp_type,
                                                       double fill_time);

/*!
 * \brief Return the velocity faces, with coefficients that add for a repeated face, in the normal-stress term at the
 * boundary face \a i normal to \a bdry_normal_axis (the lower boundary if \a bdry_is_lower). If \a bdry_interp_type is
 * QUADRATIC, the term is (11/3)*u_I - 2*u_B - (2/3)*u_I2 - u_div, less (4/3)*(u_I - u_div) if \a pseudo_traction. If it
 * is LINEAR, the term is u_I - u_div, and none if \a pseudo_traction. Here u_B, u_I, and u_I2 are the normal velocities
 * at \a i and the next two faces inside the domain, and u_div is get_divergence_free_ghost_value() in the ghost cell
 * outside the face \a i.
 */
StokesBcCoefStrategy::SideStencil get_normal_stress_stencil(const SAMRAI::hier::Index<NDIM>& i,
                                                            unsigned int bdry_normal_axis,
                                                            bool bdry_is_lower,
                                                            const double* dx,
                                                            IBTK::BdryInterpType bdry_interp_type,
                                                            bool pseudo_traction);
} // namespace detail
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
