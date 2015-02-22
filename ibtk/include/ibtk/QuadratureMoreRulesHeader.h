// this is to include all the header files on anisotropic quadrature rules in current folder: 

#ifndef included_QuadratureMoreRulesHeader
#define included_QuadratureMoreRulesHeader

// include tradition rules: newton cote

#include "ibtk/QuadratureNewtonCotes.h"
// include anisotropic rules:

#include "ibtk/QuadratureAnisotropicNewtonCotes.h"
#include "ibtk/QuadratureAnisotropicGauss.h"
#include "ibtk/QuadratureAnisotropicGrid.h"
// include anisotropic composite rules:
#include "ibtk/QuadratureAnisotropicCompositeNewtonCotes.h"
#include "ibtk/QuadratureAnisotropicCompositeGauss.h"

#endif
