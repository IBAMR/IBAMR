#ifndef included_IBAMR_RNG
#define included_IBAMR_RNG

namespace IBAMR
{
/*!
 * \brief Class RNG organizes functions that provide random-number generator
 * functionality.
 */
class RNG
{
public:
    static void srandgen(unsigned long seed);

    static void genrand(double* rn);

    static void genrandn(double* result);

    static void parallel_seed(int global_seed);

private:
    RNG() = delete;
    RNG(RNG&) = delete;
    ~RNG() = delete;
    RNG& operator=(RNG&) = delete;
};
} // namespace IBAMR

#endif //#ifndef included_IBAMR_RNG
