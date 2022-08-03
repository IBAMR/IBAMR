#ifndef included_EulerianFEStructure
#define included_EulerianFEStructure

#include <tbox/Pointer.h>

namespace IBTK
{
class FEDataManager;
}

struct EulerianFEStructure
{
    double d_R;
    IBTK::Vector3d d_X0;
    int d_chi_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_chi_var;
    IBTK::FEDataManager* d_fe_data_manager;
};

#endif // included_EulerianFEStructure