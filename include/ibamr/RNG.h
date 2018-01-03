#ifndef included_IBAMR_RNG
#define included_IBAMR_RNG

#include <petscsys.h>

namespace IBAMR
{
/*!
 * \brief Class RNG organizes functions that provide random-number generator
 * functionality.
 */
class PETSC_VISIBILITY_PUBLIC RNG
{
public:
    static void srandgen(unsigned long seed);

    static void genrand(double* rn);

    static void genrandn(double* result);

    static void parallel_seed(int global_seed);

private:
    RNG();
    RNG(RNG&);
    ~RNG();
    RNG& operator=(RNG&);
};
} // namespace IBAMR

#endif //#ifndef included_IBAMR_RNG
