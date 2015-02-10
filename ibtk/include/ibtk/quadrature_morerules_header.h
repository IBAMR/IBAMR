// this is to include all the header files in current folder: quadrature_morerules

#ifndef QUADRATRE_MORE_RULES_HEADER_H
#define QUADRATRE_MORE_RULES_HEADER_H

// include tradition rules: newton cote

#include "ibtk/quadrature_newtoncotes.h"
// include anisotropic rules:

#include "ibtk/quadrature_anisotropic_newtoncotes.h"
#include "ibtk/quadrature_anisotropic_gauss.h"
#include "ibtk/quadrature_anisotropic_grid.h"
// include anisotropic composite rules:
#include "ibtk/quadrature_anisotropic_composite_newtoncotes.h"
#include "ibtk/quadrature_anisotropic_composite_gauss.h"

#endif
