// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_PhaseChangeExamples_Applications
#define included_PhaseChangeExamples_Applications

#include <ibamr/config.h>

namespace PhaseChangeExamples
{
// Each entry owns IBTK initialization and all application objects through shutdown.
int run_stefan(int argc, char* argv[]);
int run_thermocapillary(int argc, char* argv[]);
} // namespace PhaseChangeExamples

#endif
