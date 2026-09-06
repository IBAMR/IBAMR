#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <fstream>
#include <string>

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv);

    std::ofstream out("output");

    const std::string iteration_filename = IBTK::format_iteration_output_filename(42, "output", "data");
    out << "Iteration filename:\n";
    out << iteration_filename << "\n\n";

    const std::string samrai_filename = IBTK::format_samrai_output_filename(42, "output", "data");
    out << "SAMRAI filename:\n";
    out << samrai_filename << "\n";

    return 0;
}
