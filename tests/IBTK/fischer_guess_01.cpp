// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/FischerGuess.h>

#include <libmesh/libmesh.h>
#include <libmesh/petsc_vector.h>

#include <fstream>

int
main(int argc, char** argv)
{
    using namespace libMesh;
    LibMeshInit libmesh_init(argc, argv, MPI_COMM_SELF);
    std::ofstream out("output");

    // Make sure that we project onto a single vector
    {
        IBTK::FischerGuess guess(1);
        PetscVector<double> solution(libmesh_init.comm(), 3, 3);
        PetscVector<double> rhs(libmesh_init.comm(), 3, 3);

        solution.set(0, 1);
        solution.set(1, 1);
        solution.set(2, 1);
        solution.close();

        rhs.set(0, 1);
        rhs.set(1, 1);
        rhs.close();

        guess.submit(solution, rhs);

        PetscVector<double> new_rhs(libmesh_init.comm(), 3, 3);
        new_rhs.set(0, 1);
        guess.guess(solution, new_rhs);

        solution.print(out);
    }

    // Same but use the same vector over and over again (i.e., with a singular
    // projection matrix)
    {
        IBTK::FischerGuess guess(3);

        for (int i = 0; i < 10; ++i)
        {
            PetscVector<double> solution(libmesh_init.comm(), 3, 3);
            PetscVector<double> rhs(libmesh_init.comm(), 3, 3);

            solution.set(0, 1);
            solution.set(1, 1);
            solution.set(2, 1);
            solution.close();

            rhs.set(0, 1);
            rhs.set(1, 1);
            rhs.close();

            guess.submit(solution, rhs);
        }

        PetscVector<double> new_rhs(libmesh_init.comm(), 3, 3);
        new_rhs.set(0, 1);
        new_rhs.close();
        PetscVector<double> solution(libmesh_init.comm(), 3, 3);
        guess.guess(solution, new_rhs);

        solution.print(out);
    }

    // Check with two vectors
    {
        IBTK::FischerGuess guess(10);

        PetscVector<double> solution_1(libmesh_init.comm(), 3, 3);
        PetscVector<double> rhs_1(libmesh_init.comm(), 3, 3);

        solution_1.set(0, 1);
        solution_1.set(1, 1);
        solution_1.set(2, 1);
        solution_1.close();

        rhs_1.set(0, 1);
        rhs_1.set(1, 1);
        rhs_1.close();
        guess.submit(solution_1, rhs_1);

        PetscVector<double> solution_2(libmesh_init.comm(), 3, 3);
        PetscVector<double> rhs_2(libmesh_init.comm(), 3, 3);

        solution_2.set(0, 1);
        solution_2.set(1, 2);
        solution_2.set(2, 3);
        solution_2.close();

        rhs_2.set(0, 1);
        rhs_2.set(1, -1);
        rhs_2.close();

        guess.submit(solution_2, rhs_2);

        PetscVector<double> new_rhs(libmesh_init.comm(), 3, 3);
        new_rhs.set(0, 1);
        PetscVector<double> solution(libmesh_init.comm(), 3, 3);
        guess.guess(solution, new_rhs);

        solution.print(out);
    }

    // Do the same thing but duplicate the vectors
    {
        IBTK::FischerGuess guess(10);

        for (int i = 0; i < 10; ++i)
        {
            PetscVector<double> solution_1(libmesh_init.comm(), 3, 3);
            PetscVector<double> rhs_1(libmesh_init.comm(), 3, 3);

            solution_1.set(0, 1);
            solution_1.set(1, 1);
            solution_1.set(2, 1);
            solution_1.close();

            rhs_1.set(0, 1);
            rhs_1.set(1, 1);
            rhs_1.close();
            guess.submit(solution_1, rhs_1);
        }

        for (int i = 0; i < 9; ++i)
        {
            PetscVector<double> solution_2(libmesh_init.comm(), 3, 3);
            PetscVector<double> rhs_2(libmesh_init.comm(), 3, 3);

            solution_2.set(0, 1);
            solution_2.set(1, 2);
            solution_2.set(2, 3);
            solution_2.close();

            rhs_2.set(0, 1);
            rhs_2.set(1, -1);
            rhs_2.close();

            guess.submit(solution_2, rhs_2);
        }

        PetscVector<double> new_rhs(libmesh_init.comm(), 3, 3);
        new_rhs.set(0, 1);
        PetscVector<double> solution(libmesh_init.comm(), 3, 3);
        guess.guess(solution, new_rhs);

        solution.print(out);
    }

    // Do something triangular
    {
        IBTK::FischerGuess guess(10);

        for (int i = 0; i < 5; ++i)
        {
            PetscVector<double> solution(libmesh_init.comm(), 10, 10);
            PetscVector<double> rhs(libmesh_init.comm(), 10, 10);

            for (int j = 0; j < i; ++j)
            {
                solution.set(j, j);
                rhs.set(j, j);
            }
            solution.close();
            rhs.close();

            guess.submit(solution, rhs);
        }

        PetscVector<double> new_rhs(libmesh_init.comm(), 10, 10);
        new_rhs.set(0, 1);
        new_rhs.set(1, 2);
        new_rhs.set(2, 3);
        new_rhs.set(3, 4);
        new_rhs.set(4, 5);
        PetscVector<double> solution(libmesh_init.comm(), 10, 10);
        guess.guess(solution, new_rhs);

        solution.print(out);
    }
}
