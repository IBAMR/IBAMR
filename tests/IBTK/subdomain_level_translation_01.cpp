// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/samrai_compatibility_names.h>
// SAMRAI INCLUDES
#include <ibtk/AppInitializer.h>
#include <ibtk/FEDataManager.h>
#include <ibtk/IBTKInit.h>

#include <SAMRAILogger.h>
#include <SAMRAIPointer.h>

#include <ibtk/app_namespaces.h>

#include "../tests.h"

int
main(int argc, char** argv)
{
    using namespace IBTK;
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Since this is a test we do not want to print file names or line numbers
    // to output files:
    SAMRAIPointer<SAMRAILogger::Appender> abort_append(new TestAppender());
    SAMRAILogger::getInstance()->setAbortAppender(abort_append);

    SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test.log");
    SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

    std::set<libMesh::subdomain_id_type> ids;
    for (libMesh::subdomain_id_type i = 0; i < 20; ++i) ids.insert(i);
    for (libMesh::subdomain_id_type i = 1000; i < 1010; ++i) ids.insert(i);

    SubdomainToPatchLevelTranslation trans(10, ids, input_db->getDatabase("trans"));

    for (unsigned int i = 0; i < 20; ++i)
    {
        plog << "subdomain id: " << i << " level: " << trans[i] << '\n';
    }

    for (unsigned int i = 1000; i < 1010; ++i)
    {
        plog << "subdomain id: " << i << " level: " << trans[i] << '\n';
    }
}
