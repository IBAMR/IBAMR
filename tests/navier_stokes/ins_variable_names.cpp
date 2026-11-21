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

#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>

#include <Variable.h>

#include <string>

#include <ibamr/app_namespaces.h>

namespace
{
void
print_variable_names(const std::string& label, const INSHierarchyIntegrator& integrator)
{
    plog << label << ":\n"
         << "  velocity: " << integrator.getVelocityVariable()->getName() << '\n'
         << "  pressure: " << integrator.getPressureVariable()->getName() << '\n'
         << "  body force: " << integrator.getBodyForceVariable()->getName() << '\n'
         << "  velocity divergence: " << integrator.getVelocityDivergenceVariable()->getName() << '\n';
}
} // namespace

// Print the names of the primary variables of each INS integrator, which begin with the integrator's object name.
int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app->getInputDatabase();

    INSStaggeredHierarchyIntegrator staggered(
        "StaggeredIntegrator", input_db->getDatabase("StaggeredIntegrator"), /*register_for_restart*/ false);
    print_variable_names("staggered", staggered);

    INSCollocatedHierarchyIntegrator collocated(
        "CollocatedIntegrator", input_db->getDatabase("CollocatedIntegrator"), /*register_for_restart*/ false);
    print_variable_names("collocated", collocated);

    INSVCStaggeredNonConservativeHierarchyIntegrator vc(
        "VCIntegrator", input_db->getDatabase("VCIntegrator"), /*register_for_restart*/ false);
    print_variable_names("variable coefficient", vc);
    return 0;
}
