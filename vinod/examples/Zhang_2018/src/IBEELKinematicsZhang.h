// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBEELKinematicsZhang
#define included_IBEELKinematicsZhang

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "../../IBEELKinematics.h"

#include <tbox/Database.h>
#include <tbox/Pointer.h>

namespace IBAMR
{
/*!
 * \brief Class IBEELKinematicsZhang implements strict Zhang et al. (2018) kinematics.
 *
 * This class enforces FIXED kinematics as specified in:
 * Zhang, C., Huang, H., & Lu, X.-Y. (2018). Effects of Reynolds number and
 * thickness on an undulatory self-propelled foil. Physics of Fluids, 30, 071902.
 *
 * Key differences from adaptive IBEELKinematics:
 * - NO Reynolds number dependent amplitude adaptation
 * - NO thickness dependent frequency adaptation
 * - NO swimming mode dependent envelope changes
 * - FIXED parameters: A_max = 0.125, linear envelope, Re varied only via viscosity
 *
 * This class guarantees compliance with Zhang's experimental protocol for
 * validation and comparison studies.
 */
class IBEELKinematicsZhang : public IBEELKinematics
{
public:
    /*!
     * \brief Constructor.
     */
    IBEELKinematicsZhang(const std::string& object_name,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                         IBTK::LDataManager* l_data_manager,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                         bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~IBEELKinematicsZhang();

protected:
    /*!
     * \brief Override: Calculate adaptive kinematics (DISABLED for Zhang mode).
     *
     * Zhang et al. (2018) uses FIXED kinematics across all Re and thickness.
     * This override enforces:
     *   - d_adapted_amplitude = 0.125 (constant)
     *   - d_adapted_frequency = base_frequency (constant)
     *   - d_envelope_power = 1.0 (linear, constant)
     *   - d_adapted_wavelength = 1.0 (constant)
     *
     * No Re-dependent or thickness-dependent adaptation is performed.
     */
    virtual void calculateAdaptiveKinematics(const double time) override;

private:
    /*!
     * \brief Copy constructor (not implemented).
     */
    IBEELKinematicsZhang(const IBEELKinematicsZhang& from);

    /*!
     * \brief Assignment operator (not implemented).
     */
    IBEELKinematicsZhang& operator=(const IBEELKinematicsZhang& that);

    /*!
     * \brief Verify Zhang compliance at initialization.
     */
    void verifyZhangCompliance();

    // Zhang-specific constants (enforced)
    static constexpr double ZHANG_A_MAX = 0.125;
    static constexpr double ZHANG_ENVELOPE_C0 = 0.03125;
    static constexpr double ZHANG_ENVELOPE_C1 = 1.03125;
    static constexpr double ZHANG_ENVELOPE_POWER = 1.0;
    static constexpr double ZHANG_WAVELENGTH = 1.0;

}; // IBEELKinematicsZhang

} // namespace IBAMR

#endif // #ifndef included_IBEELKinematicsZhang
