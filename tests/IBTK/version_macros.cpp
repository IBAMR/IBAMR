#include <ibtk/config.h>

#include <ibamr/config.h>

#include <fstream>

// Verify that the version check macros produce reasonable output

int main()
{
    std::ofstream out("output");

    // IBTK_VERSION_GTE is only in CMake so test for that
#ifdef IBTK_VERSION_GTE
#if IBTK_VERSION_GTE(99, 0, 0)
    out << "IBTK version greater than 99.0.0\n";
#else
    out << "IBTK version less than 99.0.0\n";
#endif

#if IBTK_VERSION_GTE(0, 0, 0)
    out << "IBTK version greater than 0.0.0\n";
#else
    out << "IBTK version less than 0.0.0\n";
#endif
#else
    // just print the right answers
    out << "IBTK version less than 99.0.0\n";
    out << "IBTK version greater than 0.0.0\n";
#endif

#if IBAMR_VERSION_GTE(99, 0, 0)
    out << "IBAMR version greater than 99.0.0\n";
#else
    out << "IBAMR version less than 99.0.0\n";
#endif

#if IBAMR_VERSION_GTE(0, 0, 0)
    out << "IBAMR version greater than 0.0.0\n";
#else
    out << "IBAMR version less than 0.0.0\n";
#endif
}
