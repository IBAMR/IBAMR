// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <fstream>

#include "cav_raw_operator_comparator.h"

int
main(int argc, char* argv[])
{
    if (argc < 2) return 1;
    std::ifstream input(argv[1]);
    if (!input) return 1;

    const IBAMR::TestSupport::CAVRawComparatorControlResults results =
        IBAMR::TestSupport::runCAVRawComparatorControls();
    const auto write_bool = [](std::ofstream& output, const char* name, const bool value)
    { output << name << " = " << (value ? "true" : "false") << "\n"; };
    std::ofstream output("output");
    write_bool(output, "exact_control", results.exact_control);
    write_bool(output, "declared_permutation_sign_gauge", results.declared_permutation_sign_gauge);
    write_bool(output, "undeclared_permutation_detected", results.undeclared_permutation_detected);
    write_bool(output, "velocity_sign_detected", results.velocity_sign_detected);
    write_bool(output, "pressure_sign_detected", results.pressure_sign_detected);
    write_bool(output, "legal_pressure_gauge", results.legal_pressure_gauge);
    write_bool(output, "illegal_velocity_shift_detected", results.illegal_velocity_shift_detected);
    write_bool(output, "omitted_dof_detected", results.omitted_dof_detected);
    write_bool(output, "omitted_row_detected", results.omitted_row_detected);
    write_bool(output, "ordered_omission_detected", results.ordered_omission_detected);
    write_bool(output, "ordered_reordering_detected", results.ordered_reordering_detected);

    const int failures =
        !(results.exact_control && results.declared_permutation_sign_gauge && results.undeclared_permutation_detected &&
          results.velocity_sign_detected && results.pressure_sign_detected && results.legal_pressure_gauge &&
          results.illegal_velocity_shift_detected && results.omitted_dof_detected && results.omitted_row_detected &&
          results.ordered_omission_detected && results.ordered_reordering_detected);
    output << "test_failures = " << failures << "\n";
    return failures;
}
