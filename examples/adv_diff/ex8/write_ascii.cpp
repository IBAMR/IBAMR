#include <mpi.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
int
main(int argc, char** argv)
{
    // setting MPI environment
    int rank, nprocs, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const char* filename = argv[1];
    MPI_File fh;
    MPI_Status status;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    MPI_File_read(fh, &size, 1, MPI_INT, &status);
    MPI_File_seek(fh, sizeof(int), MPI_SEEK_SET);
    std::vector<double> a(size);
    MPI_File_read(fh, &a[0], size, MPI_DOUBLE, &status);
    if (rank == 0)
    {
        std::string ascii_filename = "temperature_profile";
        std::ofstream fout(ascii_filename.c_str(), std::ios::out);
        fout.precision(7);
        for (int i = 0; i < size;)
        {
            fout << a[i] << "\t" << a[i + 1] << "\n";
            i = i + 2;
        }
        fout.close();
    }
    MPI_Finalize();
    return 0;
}
