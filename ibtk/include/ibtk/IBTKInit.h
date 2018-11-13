// Filename: IBTKInit.h
//
// Copyright (c) 2002-2017, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBTKInit
#define included_IBTKInit

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "SAMRAI_config.h"
#include "tbox/SAMRAIManager.h"
#include "tbox/SAMRAI_MPI.h"
#include <IBTK_config.h>

#include <petscsys.h>
#ifdef IBTK_HAVE_LIBMESH
#include "libmesh/libmesh.h"
#include "libmesh/reference_counter.h"
#endif

#include "ibtk/IBTK_MPI.h"

namespace IBTK
{
/**
 * @brief Initialization for IBAMR programs.
 *
 * The IBTKInit class handles the initializations for PETSc, LibMesh, and SAMRAI. This object should be initialized at
 * the start of the main() function. The destruction of the object correctly closes the libraries.
 *
 */

class IBTKInit
{
public:
    /**
     * Constructor for IBTKInit. Initializes libraries and sets the SAMRAI world communicator.
     */
    IBTKInit(int argc,
              char** argvmake,
              IBTK_MPI::comm communicator = PETSC_COMM_WORLD,
              char* petsc_file = nullptr,
              char* petsc_help = nullptr);

    /**
     * Destructor. Closes libraries appropriately.
     */
    ~IBTKInit();

    /**
     * Get libMesh initialization object.
     */
    std::shared_ptr<libMesh::LibMeshInit> getLibMeshInit()
    {
#ifdef IBTK_HAVE_LIBMESH
        return d_libmesh_init;
#else
        TBOX_ERROR("IBTKInit::getLibMeshInit() \n"
                   << "IBAMR not compiled with libMesh!\n");
#endif
    }

private:
#ifdef IBTK_HAVE_LIBMESH
    std::shared_ptr<libMesh::LibMeshInit> d_libmesh_init;
#endif
};

} // namespace IBTK

#endif
