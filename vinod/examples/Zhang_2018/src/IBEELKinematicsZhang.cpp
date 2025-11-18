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

//////////////////////////// INCLUDES /////////////////////////////////////////

#include "IBEELKinematicsZhang.h"

#include "ibtk/IBTK_MPI.h"
#include "tbox/MathUtilities.h"

#include <iostream>
#include <iomanip>
#include <cmath>

#include "ibamr/namespaces.h"

namespace IBAMR
{

IBEELKinematicsZhang::IBEELKinematicsZhang(const std::string& object_name,
                                           Pointer<Database> input_db,
                                           LDataManager* l_data_manager,
                                           Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                           bool register_for_restart)
    : IBEELKinematics(object_name, input_db, l_data_manager, patch_hierarchy, register_for_restart)
{
    // Verify Zhang compliance and warn about any parameter mismatches
    verifyZhangCompliance();

    // Override base class defaults to enforce Zhang parameters
    d_base_amplitude = ZHANG_A_MAX;
    d_envelope_power = ZHANG_ENVELOPE_POWER;
    d_adapted_wavelength = ZHANG_WAVELENGTH;

    // Force disable adaptation (ignore input file setting)
    d_enable_shape_adaptation = false;

    if (IBTK_MPI::getRank() == 0)
    {
        std::cout << "\n======================================================" << std::endl;
        std::cout << "  Zhang et al. (2018) Strict Kinematics Mode ACTIVE" << std::endl;
        std::cout << "======================================================" << std::endl;
        std::cout << "  Reference: Physics of Fluids 30, 071902 (2018)" << std::endl;
        std::cout << "  FIXED kinematics (no Re/thickness adaptation):" << std::endl;
        std::cout << "    A_max     = " << ZHANG_A_MAX << " (constant)" << std::endl;
        std::cout << "    Envelope  = " << ZHANG_A_MAX << " × (X + " << ZHANG_ENVELOPE_C0
                  << ") / " << ZHANG_ENVELOPE_C1 << std::endl;
        std::cout << "    Power     = " << ZHANG_ENVELOPE_POWER << " (linear)" << std::endl;
        std::cout << "    Wavelength= " << ZHANG_WAVELENGTH << " (one wavelength)" << std::endl;
        std::cout << "  Current simulation parameters:" << std::endl;
        std::cout << "    Re        = " << d_reynolds_number << std::endl;
        std::cout << "    h/L       = " << d_thickness_ratio << std::endl;
        std::cout << "    Frequency = " << d_base_frequency << std::endl;
        std::cout << "======================================================\n" << std::endl;
    }

    return;

} // IBEELKinematicsZhang constructor

IBEELKinematicsZhang::~IBEELKinematicsZhang()
{
    // Destructor (base class handles cleanup)
    return;

} // ~IBEELKinematicsZhang

void
IBEELKinematicsZhang::calculateAdaptiveKinematics(const double time)
{
    // ZHANG MODE: NO ADAPTATION
    //
    // Zhang et al. (2018) uses FIXED kinematics across all Re and thickness values.
    // This is critical for their experimental protocol: they isolate the effect of
    // flow regime (Re) and geometry (thickness) by keeping the swimming pattern constant.
    //
    // Any modification of amplitude, frequency, or envelope shape based on Re or
    // thickness would violate Zhang's approach and make comparison impossible.

    // Enforce fixed Zhang parameters (no Re or thickness dependence)
    d_adapted_amplitude = ZHANG_A_MAX;           // Always 0.125
    d_adapted_frequency = d_base_frequency;      // Use base frequency (no adaptation)
    d_adapted_wavelength = ZHANG_WAVELENGTH;     // Always 1.0
    d_envelope_power = ZHANG_ENVELOPE_POWER;     // Always 1.0 (linear envelope)

    // Log parameters periodically (but less frequently than adaptive mode)
    static bool first_call = true;
    const double log_interval = 2.0;  // Log every 2 time units
    static double last_log_time = -log_interval;

    if (first_call || (time - last_log_time) >= log_interval)
    {
        if (IBTK_MPI::getRank() == 0)
        {
            std::cout << "\n=== Zhang (2018) Kinematics Check (t=" << time << ") ===" << std::endl;
            std::cout << "  FIXED parameters (no adaptation):" << std::endl;
            std::cout << "    Amplitude   = " << d_adapted_amplitude << " (constant)" << std::endl;
            std::cout << "    Frequency   = " << d_adapted_frequency << " (constant)" << std::endl;
            std::cout << "    Envelope    = linear (power = " << d_envelope_power << ")" << std::endl;
            std::cout << "    Wavelength  = " << d_adapted_wavelength << std::endl;
            std::cout << "  Re = " << d_reynolds_number << ", h/L = " << d_thickness_ratio << std::endl;
            std::cout << "====================================================\n" << std::endl;
        }
        first_call = false;
        last_log_time = time;
    }

    return;

} // calculateAdaptiveKinematics

void
IBEELKinematicsZhang::verifyZhangCompliance()
{
    // Check if input parameters match Zhang specification
    // Warn if there are discrepancies (but don't error - we'll override them)

    bool compliance_warning = false;

    if (IBTK_MPI::getRank() == 0)
    {
        std::cout << "\n=== Zhang (2018) Compliance Verification ===" << std::endl;

        // Check amplitude
        if (std::abs(d_base_amplitude - ZHANG_A_MAX) > 1e-6)
        {
            std::cout << "  WARNING: base_amplitude = " << d_base_amplitude
                      << " (Zhang requires " << ZHANG_A_MAX << ")" << std::endl;
            std::cout << "           Will be overridden to " << ZHANG_A_MAX << std::endl;
            compliance_warning = true;
        }
        else
        {
            std::cout << "  ✓ Amplitude = " << ZHANG_A_MAX << " (correct)" << std::endl;
        }

        // Check if adaptation was enabled in input file
        if (d_enable_shape_adaptation)
        {
            std::cout << "  WARNING: enable_shape_adaptation = TRUE in input file" << std::endl;
            std::cout << "           Will be forced to FALSE for Zhang mode" << std::endl;
            compliance_warning = true;
        }
        else
        {
            std::cout << "  ✓ Shape adaptation disabled (correct)" << std::endl;
        }

        // Check envelope power
        if (std::abs(d_envelope_power - ZHANG_ENVELOPE_POWER) > 1e-6)
        {
            std::cout << "  WARNING: envelope_power = " << d_envelope_power
                      << " (Zhang requires " << ZHANG_ENVELOPE_POWER << ")" << std::endl;
            std::cout << "           Will be overridden to " << ZHANG_ENVELOPE_POWER << std::endl;
            compliance_warning = true;
        }
        else
        {
            std::cout << "  ✓ Envelope power = " << ZHANG_ENVELOPE_POWER << " (linear, correct)" << std::endl;
        }

        // Check swimming mode (should be 0.0 for anguilliform)
        if (std::abs(d_swimming_mode) > 1e-6)
        {
            std::cout << "  WARNING: swimming_mode = " << d_swimming_mode
                      << " (Zhang uses 0.0 for anguilliform)" << std::endl;
            compliance_warning = true;
        }
        else
        {
            std::cout << "  ✓ Swimming mode = 0.0 (anguilliform, correct)" << std::endl;
        }

        // Verify Re is within Zhang's test range
        const double zhang_re_min = 50.0;
        const double zhang_re_max = 200000.0;
        if (d_reynolds_number < zhang_re_min || d_reynolds_number > zhang_re_max)
        {
            std::cout << "  WARNING: Re = " << d_reynolds_number
                      << " is outside Zhang's range [" << zhang_re_min
                      << ", " << zhang_re_max << "]" << std::endl;
            compliance_warning = true;
        }
        else
        {
            std::cout << "  ✓ Re = " << d_reynolds_number << " (within Zhang range)" << std::endl;
        }

        // Verify thickness is within Zhang's test range
        const double zhang_h_min = 0.04;
        const double zhang_h_max = 0.24;
        if (d_thickness_ratio < zhang_h_min || d_thickness_ratio > zhang_h_max)
        {
            std::cout << "  WARNING: h/L = " << d_thickness_ratio
                      << " is outside Zhang's range [" << zhang_h_min
                      << ", " << zhang_h_max << "]" << std::endl;
            compliance_warning = true;
        }
        else
        {
            std::cout << "  ✓ h/L = " << d_thickness_ratio << " (within Zhang range)" << std::endl;
        }

        if (compliance_warning)
        {
            std::cout << "\n  NOTE: Some parameters differ from Zhang specification." << std::endl;
            std::cout << "        Critical parameters will be overridden automatically." << std::endl;
            std::cout << "        For exact Zhang reproduction, update your input file.\n" << std::endl;
        }
        else
        {
            std::cout << "\n  ✓ All parameters comply with Zhang (2018) specification.\n" << std::endl;
        }

        std::cout << "============================================\n" << std::endl;
    }

    return;

} // verifyZhangCompliance

} // namespace IBAMR
