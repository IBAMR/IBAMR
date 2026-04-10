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

#ifndef included_IBAMR_INSTurbulenceModel
#define included_IBAMR_INSTurbulenceModel

#include <ibamr/config.h>

#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

#include <string>
#include <utility>
#include <vector>

namespace IBAMR
{
class StokesSpecifications;
}
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class SideVariable;
} // namespace pdat
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

namespace IBAMR
{
/*!
 * \brief Abstract interface for turbulence closures coupled to
 * INSStaggeredHierarchyIntegrator.
 *
 * Implementations are responsible for computing an explicit side-centered
 * forcing contribution to the momentum equation from a velocity field with
 * valid ghost-cell values.
 */
class INSTurbulenceModel : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    explicit INSTurbulenceModel(std::string object_name) : d_object_name(std::move(object_name))
    {
    }

    /*!
     * \brief Destructor.
     */
    ~INSTurbulenceModel() override = default;

    /*!
     * Compute the turbulence-model contribution to the momentum equation.
     *
     * Implementations should overwrite \p F_idx on the supplied hierarchy.
     * The velocity data in \p U_idx are assumed to already have valid ghost
     * cells.
     */
    virtual void computeTurbulenceForce(int F_idx,
                                        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> F_var,
                                        int U_idx,
                                        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> U_var,
                                        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                        double data_time,
                                        const StokesSpecifications& problem_coefs) = 0;

protected:
    /*!
     * Object name used in logging and subordinate database/variable names.
     */
    std::string d_object_name;
};
} // namespace IBAMR

#endif
