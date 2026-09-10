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

#ifndef included_IBAMR_INSTurbulenceStatistics
#define included_IBAMR_INSTurbulenceStatistics

#include <ibamr/config.h>

#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

#include <string>
#include <utility>
#include <vector>

namespace IBTK
{
class HierarchyMathOps;
}
namespace SAMRAI
{
namespace appu
{
template <int DIM>
class VisItDataWriter;
} // namespace appu
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
 * \brief Abstract interface for turbulence-statistics accumulation coupled to
 * INSStaggeredHierarchyIntegrator.
 */
class INSTurbulenceStatistics : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name Name used for diagnostics and subordinate objects.
     * \param statistics_start_time Simulation time at which statistics
     * accumulation should begin.
     */
    explicit INSTurbulenceStatistics(std::string object_name, double statistics_start_time = 0.0);

    /*!
     * \brief Destructor.
     */
    ~INSTurbulenceStatistics() override = default;

    /*!
     * Update stored statistics from the supplied velocity field.
     *
     * Returns true when the tracked statistics report that they are at
     * periodic steady state.
     */
    virtual bool updateStatistics(int U_idx,
                                  SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> U_var,
                                  const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                                  double data_time,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                  SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops) = 0;

    /*!
     * \brief Return whether the tracked statistics are at periodic steady
     * state.
     */
    virtual bool isAtSteadyState() const = 0;

    /*!
     * \brief Register the VisIt data writer used by the hierarchy
     * integrator.
     *
     * The default implementation is a no-op.
     */
    virtual void registerVisItDataWriter(SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM>> visit_writer);

    /*!
     * \brief Prepare any statistic-specific plot data on the supplied
     * hierarchy at the specified output time.
     *
     * The default implementation is a no-op.
     */
    virtual void setupPlotData(double data_time,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops);

    /*!
     * \brief Deallocate any statistic-specific plot data owned by this object.
     *
     * The default implementation is a no-op.
     */
    virtual void deallocatePlotData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

protected:
    /*!
     * \brief Return true when statistics should be updated at the supplied
     * simulation time.
     */
    virtual bool shouldUpdateStatistics(double data_time);

    /*!
     * Object name used in logging and subordinate database/variable names.
     */
    std::string d_object_name;

    /*!
     * Simulation time at which statistics accumulation should begin.
     */
    double d_statistics_start_time = 0.0;
};
} // namespace IBAMR

#include <ibamr/private/INSTurbulenceStatistics-inl.h> // IWYU pragma: keep

#endif
