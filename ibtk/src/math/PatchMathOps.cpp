// Filename: PatchMathOps.cpp
// Created on 23 Jul 2002 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ostream>

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "EdgeData.h" // IWYU pragma: keep
#include "FaceData.h"
#include "FaceGeometry.h"
#include "IBTK_config.h"
#include "IntVector.h"
#include "NodeData.h"
#include "NodeGeometry.h"
#include "Patch.h"
#include "PatchFaceDataOpsReal.h"
#include "PatchSideDataOpsReal.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "boost/array.hpp"
#include "ibtk/PatchMathOps.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define LAPLACE_FC IBTK_FC_FUNC(laplace2d, LAPLACE2D)
#define LAPLACE_ADD_FC IBTK_FC_FUNC(laplaceadd2d, LAPLACEADD2D)
#define DAMPED_LAPLACE_FC IBTK_FC_FUNC(dampedlaplace2d, DAMPEDLAPLACE2D)
#define DAMPED_LAPLACE_ADD_FC IBTK_FC_FUNC(dampedlaplaceadd2d, DAMPEDLAPLACEADD2D)

#define MULTIPLY1_FC IBTK_FC_FUNC(multiply12d, MULTIPLY12D)
#define MULTIPLY_ADD1_FC IBTK_FC_FUNC(multiplyadd12d, MULTIPLYADD12D)
#define MULTIPLY2_FC IBTK_FC_FUNC(multiply22d, MULTIPLY22D)
#define MULTIPLY_ADD2_FC IBTK_FC_FUNC(multiplyadd22d, MULTIPLYADD22D)
#define MULTIPLY_ADD3_FC IBTK_FC_FUNC(multiplyadd32d, MULTIPLYADD32D)

#define PW_L1_NORM_FC IBTK_FC_FUNC(pwl1norm2d, PWL1NORM2D)
#define PW_L2_NORM_FC IBTK_FC_FUNC(pwl2norm2d, PWL2NORM2D)
#define PW_MAX_NORM_FC IBTK_FC_FUNC(pwmaxnorm2d, PWMAXNORM2D)

#define C_TO_C_CURL_FC IBTK_FC_FUNC(ctoccurl2d, CTOCCURL2D)
#define C_TO_C_DIV_FC IBTK_FC_FUNC(ctocdiv2d, CTOCDIV2D)
#define C_TO_C_DIV_ADD_FC IBTK_FC_FUNC(ctocdivadd2d, CTOCDIVADD2D)
#define C_TO_C_GRAD_FC IBTK_FC_FUNC(ctocgrad2d, CTOCGRAD2D)
#define C_TO_C_GRAD_ADD_FC IBTK_FC_FUNC(ctocgradadd2d, CTOCGRADADD2D)

#define C_TO_C_ANISO_F_LAPLACE_FC IBTK_FC_FUNC(ctocanisoflaplace2d, CTOCANISOFLAPLACE2D)
#define C_TO_C_ANISO_F_LAPLACE_ADD_FC IBTK_FC_FUNC(ctocanisoflaplaceadd2d, CTOCANISOFLAPLACEADD2D)
#define C_TO_C_ANISO_F_DAMPED_LAPLACE_FC IBTK_FC_FUNC(ctocanisofdampedlaplace2d, CTOCANISOFDAMPEDLAPLACE2D)
#define C_TO_C_ANISO_F_DAMPED_LAPLACE_ADD_FC IBTK_FC_FUNC(ctocanisofdampedlaplaceadd2d, CTOCANISOFDAMPEDLAPLACEADD2D)

#define C_TO_C_ANISO_S_LAPLACE_FC IBTK_FC_FUNC(ctocanisoslaplace2d, CTOCANISOSLAPLACE2D)
#define C_TO_C_ANISO_S_LAPLACE_ADD_FC IBTK_FC_FUNC(ctocanisoslaplaceadd2d, CTOCANISOSLAPLACEADD2D)
#define C_TO_C_ANISO_S_DAMPED_LAPLACE_FC IBTK_FC_FUNC(ctocanisosdampedlaplace2d, CTOCANISOSDAMPEDLAPLACE2D)
#define C_TO_C_ANISO_S_DAMPED_LAPLACE_ADD_FC IBTK_FC_FUNC(ctocanisosdampedlaplaceadd2d, CTOCANISOSDAMPEDLAPLACEADD2D)

#define C_TO_F_GRAD_FC IBTK_FC_FUNC(ctofgrad2d, CTOFGRAD2D)
#define C_TO_F_FLUX_FC IBTK_FC_FUNC(ctofflux2d, CTOFFLUX2D)
#define C_TO_F_ANISO_FLUX_FC IBTK_FC_FUNC(ctofanisoflux2d, CTOFANISOFLUX2D)
#define C_TO_F_GRAD_ADD_FC IBTK_FC_FUNC(ctofgradadd2d, CTOFGRADADD2D)
#define C_TO_F_INTERP_FC IBTK_FC_FUNC(ctofinterp2nd2d, CTOFINTERP2ND2D)

#define C_TO_S_GRAD_FC IBTK_FC_FUNC(ctosgrad2d, CTOSGRAD2D)
#define C_TO_S_FLUX_FC IBTK_FC_FUNC(ctosflux2d, CTOSFLUX2D)
#define C_TO_S_ANISO_FLUX_FC IBTK_FC_FUNC(ctosanisoflux2d, CTOSANISOFLUX2D)
#define C_TO_S_GRAD_ADD_FC IBTK_FC_FUNC(ctosgradadd2d, CTOSGRADADD2D)
#define C_TO_S_INTERP_FC IBTK_FC_FUNC(ctosinterp2nd2d, CTOSINTERP2ND2D)

#define F_TO_C_CURL_FC IBTK_FC_FUNC(ftoccurl2d, FTOCCURL2D)
#define F_TO_C_DIV_FC IBTK_FC_FUNC(ftocdiv2d, FTOCDIV2D)
#define F_TO_C_DIV_ADD_FC IBTK_FC_FUNC(ftocdivadd2d, FTOCDIVADD2D)
#define F_TO_C_INTERP_FC IBTK_FC_FUNC(ftocinterp2nd2d, FTOCINTERP2ND2D)

#define S_TO_C_CURL_FC IBTK_FC_FUNC(stoccurl2d, STOCCURL2D)
#define S_TO_C_DIV_FC IBTK_FC_FUNC(stocdiv2d, STOCDIV2D)
#define S_TO_C_DIV_ADD_FC IBTK_FC_FUNC(stocdivadd2d, STOCDIVADD2D)
#define S_TO_C_INTERP_FC IBTK_FC_FUNC(stocinterp2nd2d, STOCINTERP2ND2D)

#define S_TO_S_VC_LAPLACE_FC IBTK_FC_FUNC(stosvclaplace2d, STOSVCLAPLACE2D)

#define N_TO_S_ROT_FC IBTK_FC_FUNC(ntosrot2d, NTOSROT2D)
#define C_TO_S_ROT_FC IBTK_FC_FUNC(ctosrot2d, CTOSROT2D)

#define S_TO_N_CURL_FC IBTK_FC_FUNC(stoncurl2d, STONCURL2D)
#endif // if (NDIM == 2)

#if (NDIM == 3)
#define LAPLACE_FC IBTK_FC_FUNC(laplace3d, LAPLACE3D)
#define LAPLACE_ADD_FC IBTK_FC_FUNC(laplaceadd3d, LAPLACEADD3D)
#define DAMPED_LAPLACE_FC IBTK_FC_FUNC(dampedlaplace3d, DAMPEDLAPLACE3D)
#define DAMPED_LAPLACE_ADD_FC IBTK_FC_FUNC(dampedlaplaceadd3d, DAMPEDLAPLACEADD3D)

#define MULTIPLY1_FC IBTK_FC_FUNC(multiply13d, MULTIPLY13D)
#define MULTIPLY_ADD1_FC IBTK_FC_FUNC(multiplyadd13d, MULTIPLYADD13D)
#define MULTIPLY2_FC IBTK_FC_FUNC(multiply23d, MULTIPLY23D)
#define MULTIPLY_ADD2_FC IBTK_FC_FUNC(multiplyadd23d, MULTIPLYADD23D)
#define MULTIPLY_ADD3_FC IBTK_FC_FUNC(multiplyadd33d, MULTIPLYADD33D)

#define PW_L1_NORM_FC IBTK_FC_FUNC(pwl1norm3d, PWL1NORM3D)
#define PW_L2_NORM_FC IBTK_FC_FUNC(pwl2norm3d, PWL2NORM3D)
#define PW_MAX_NORM_FC IBTK_FC_FUNC(pwmaxnorm3d, PWMAXNORM3D)

#define C_TO_C_CURL_FC IBTK_FC_FUNC(ctoccurl3d, CTOCCURL3D)
#define C_TO_C_DIV_FC IBTK_FC_FUNC(ctocdiv3d, CTOCDIV3D)
#define C_TO_C_DIV_ADD_FC IBTK_FC_FUNC(ctocdivadd3d, CTOCDIVADD3D)
#define C_TO_C_GRAD_FC IBTK_FC_FUNC(ctocgrad3d, CTOCGRAD3D)
#define C_TO_C_GRAD_ADD_FC IBTK_FC_FUNC(ctocgradadd3d, CTOCGRADADD3D)

#define C_TO_C_ANISO_F_LAPLACE_FC IBTK_FC_FUNC(ctocanisoflaplace3d, CTOCANISOFLAPLACE3D)
#define C_TO_C_ANISO_F_LAPLACE_ADD_FC IBTK_FC_FUNC(ctocanisoflaplaceadd3d, CTOCANISOFLAPLACEADD3D)
#define C_TO_C_ANISO_F_DAMPED_LAPLACE_FC IBTK_FC_FUNC(ctocanisofdampedlaplace3d, CTOCANISOFDAMPEDLAPLACE3D)
#define C_TO_C_ANISO_F_DAMPED_LAPLACE_ADD_FC IBTK_FC_FUNC(ctocanisofdampedlaplaceadd3d, CTOCANISOFDAMPEDLAPLACEADD3D)

#define C_TO_C_ANISO_S_LAPLACE_FC IBTK_FC_FUNC(ctocanisoslaplace3d, CTOCANISOSLAPLACE3D)
#define C_TO_C_ANISO_S_LAPLACE_ADD_FC IBTK_FC_FUNC(ctocanisoslaplaceadd3d, CTOCANISOSLAPLACEADD3D)
#define C_TO_C_ANISO_S_DAMPED_LAPLACE_FC IBTK_FC_FUNC(ctocanisosdampedlaplace3d, CTOCANISOSDAMPEDLAPLACE3D)
#define C_TO_C_ANISO_S_DAMPED_LAPLACE_ADD_FC IBTK_FC_FUNC(ctocanisosdampedlaplaceadd3d, CTOCANISOSDAMPEDLAPLACEADD3D)

#define C_TO_F_GRAD_FC IBTK_FC_FUNC(ctofgrad3d, CTOFGRAD3D)
#define C_TO_F_FLUX_FC IBTK_FC_FUNC(ctofflux3d, CTOFFLUX3D)
#define C_TO_F_ANISO_FLUX_FC IBTK_FC_FUNC(ctofanisoflux3d, CTOFANISOFLUX3D)
#define C_TO_F_GRAD_ADD_FC IBTK_FC_FUNC(ctofgradadd3d, CTOFGRADADD3D)
#define C_TO_F_INTERP_FC IBTK_FC_FUNC(ctofinterp2nd3d, CTOFINTERP2ND3D)

#define C_TO_S_GRAD_FC IBTK_FC_FUNC(ctosgrad3d, CTOSGRAD3D)
#define C_TO_S_FLUX_FC IBTK_FC_FUNC(ctosflux3d, CTOSFLUX3D)
#define C_TO_S_ANISO_FLUX_FC IBTK_FC_FUNC(ctosanisoflux3d, CTOSANISOFLUX3D)
#define C_TO_S_GRAD_ADD_FC IBTK_FC_FUNC(ctosgradadd3d, CTOSGRADADD3D)
#define C_TO_S_INTERP_FC IBTK_FC_FUNC(ctosinterp2nd3d, CTOSINTERP2ND3D)

#define F_TO_C_CURL_FC IBTK_FC_FUNC(ftoccurl3d, FTOCCURL3D)
#define F_TO_C_DIV_FC IBTK_FC_FUNC(ftocdiv3d, FTOCDIV3D)
#define F_TO_C_DIV_ADD_FC IBTK_FC_FUNC(ftocdivadd3d, FTOCDIVADD3D)
#define F_TO_C_INTERP_FC IBTK_FC_FUNC(ftocinterp2nd3d, FTOCINTERP2ND3D)

#define F_TO_F_CURL_FC IBTK_FC_FUNC(ftofcurl3d, FTOFCURL3D)

#define S_TO_C_CURL_FC IBTK_FC_FUNC(stoccurl3d, STOCCURL3D)
#define S_TO_C_DIV_FC IBTK_FC_FUNC(stocdiv3d, STOCDIV3D)
#define S_TO_C_DIV_ADD_FC IBTK_FC_FUNC(stocdivadd3d, STOCDIVADD3D)
#define S_TO_C_INTERP_FC IBTK_FC_FUNC(stocinterp2nd3d, STOCINTERP2ND3D)

#define S_TO_S_CURL_FC IBTK_FC_FUNC(stoscurl3d, STOSCURL3D)

#define S_TO_E_CURL_FC IBTK_FC_FUNC(stoecurl3d, STOECURL3D)

#define E_TO_S_ROT_FC IBTK_FC_FUNC(etosrot3d, ETOSROT3D)
#endif // if (NDIM == 3)

extern "C" {
void LAPLACE_FC(double* F,
                const int& F_gcw,
                const double& alpha,
                const double* U,
                const int& U_gcw,
                const int& ilower0,
                const int& iupper0,
                const int& ilower1,
                const int& iupper1,
#if (NDIM == 3)
                const int& ilower2,
                const int& iupper2,
#endif
                const double* dx);

void LAPLACE_ADD_FC(double* F,
                    const int& F_gcw,
                    const double& alpha,
                    const double* U,
                    const int& U_gcw,
                    const double& beta,
                    const double* V,
                    const int& V_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
#if (NDIM == 3)
                    const int& ilower2,
                    const int& iupper2,
#endif
                    const double* dx);

void DAMPED_LAPLACE_FC(double* F,
                       const int& F_gcw,
                       const double& alpha,
                       const double& beta,
                       const double* U,
                       const int& U_gcw,
                       const int& ilower0,
                       const int& iupper0,
                       const int& ilower1,
                       const int& iupper1,
#if (NDIM == 3)
                       const int& ilower2,
                       const int& iupper2,
#endif
                       const double* dx);

void DAMPED_LAPLACE_ADD_FC(double* F,
                           const int& F_gcw,
                           const double& alpha,
                           const double& beta,
                           const double* U,
                           const int& U_gcw,
                           const double& gamma,
                           const double* V,
                           const int& V_gcw,
                           const int& ilower0,
                           const int& iupper0,
                           const int& ilower1,
                           const int& iupper1,
#if (NDIM == 3)
                           const int& ilower2,
                           const int& iupper2,
#endif
                           const double* dx);

void C_TO_C_CURL_FC(double* W,
                    const int& W_gcw,
                    const double* U,
                    const int& U_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
#if (NDIM == 3)
                    const int& ilower2,
                    const int& iupper2,
#endif
                    const double* dx);

void C_TO_C_DIV_FC(double* D,
                   const int& D_gcw,
                   const double& alpha,
                   const double* U,
                   const int& U_gcw,
                   const int& ilower0,
                   const int& iupper0,
                   const int& ilower1,
                   const int& iupper1,
#if (NDIM == 3)
                   const int& ilower2,
                   const int& iupper2,
#endif
                   const double* dx);

void C_TO_C_DIV_ADD_FC(double* D,
                       const int& D_gcw,
                       const double& alpha,
                       const double* U,
                       const int& U_gcw,
                       const double& beta,
                       const double* V,
                       const int& V_gcw,
                       const int& ilower0,
                       const int& iupper0,
                       const int& ilower1,
                       const int& iupper1,
#if (NDIM == 3)
                       const int& ilower2,
                       const int& iupper2,
#endif
                       const double* dx);

void C_TO_C_GRAD_FC(double* G,
                    const int& G_gcw,
                    const double& alpha,
                    const double* U,
                    const int& U_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
#if (NDIM == 3)
                    const int& ilower2,
                    const int& iupper2,
#endif
                    const double* dx);

void C_TO_C_GRAD_ADD_FC(double* G,
                        const int& G_gcw,
                        const double& alpha,
                        const double* U,
                        const int& U_gcw,
                        const double& beta,
                        const double* V,
                        const int& V_gcw,
                        const int& ilower0,
                        const int& iupper0,
                        const int& ilower1,
                        const int& iupper1,
#if (NDIM == 3)
                        const int& ilower2,
                        const int& iupper2,
#endif
                        const double* dx);

void MULTIPLY1_FC(double* U,
                  const int& U_gcw,
                  const double& alpha,
                  const double* V,
                  const int& V_gcw,
                  const int& ilower0,
                  const int& iupper0,
                  const int& ilower1,
                  const int& iupper1
#if (NDIM == 3)
                  ,
                  const int& ilower2,
                  const int& iupper2
#endif
                  );

void MULTIPLY_ADD1_FC(double* U,
                      const int& U_gcw,
                      const double& alpha,
                      const double* V,
                      const int& V_gcw,
                      const double& beta,
                      const double* W,
                      const int& W_gcw,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1
#if (NDIM == 3)
                      ,
                      const int& ilower2,
                      const int& iupper2
#endif
                      );

void MULTIPLY2_FC(double* U,
                  const int& U_gcw,
                  const double* A,
                  const int& A_gcw,
                  const double* V,
                  const int& V_gcw,
                  const int& ilower0,
                  const int& iupper0,
                  const int& ilower1,
                  const int& iupper1
#if (NDIM == 3)
                  ,
                  const int& ilower2,
                  const int& iupper2
#endif
                  );

void MULTIPLY_ADD2_FC(double* U,
                      const int& U_gcw,
                      const double* A,
                      const int& A_gcw,
                      const double* V,
                      const int& V_gcw,
                      const double& beta,
                      const double* W,
                      const int& W_gcw,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1
#if (NDIM == 3)
                      ,
                      const int& ilower2,
                      const int& iupper2
#endif
                      );

void MULTIPLY_ADD3_FC(double* U,
                      const int& U_gcw,
                      const double* A,
                      const int& A_gcw,
                      const double* V,
                      const int& V_gcw,
                      const double* B,
                      const int& B_gcw,
                      const double* W,
                      const int& W_gcw,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1
#if (NDIM == 3)
                      ,
                      const int& ilower2,
                      const int& iupper2
#endif
                      );

void PW_L1_NORM_FC(double* U,
                   const int& U_gcw,
                   const double* V,
                   const int& V_gcw,
                   const int& V_depth,
                   const int& ilower0,
                   const int& iupper0,
                   const int& ilower1,
                   const int& iupper1
#if (NDIM == 3)
                   ,
                   const int& ilower2,
                   const int& iupper2
#endif
                   );

void PW_L2_NORM_FC(double* U,
                   const int& U_gcw,
                   const double* V,
                   const int& V_gcw,
                   const int& V_depth,
                   const int& ilower0,
                   const int& iupper0,
                   const int& ilower1,
                   const int& iupper1
#if (NDIM == 3)
                   ,
                   const int& ilower2,
                   const int& iupper2
#endif
                   );

void PW_MAX_NORM_FC(double* U,
                    const int& U_gcw,
                    const double* V,
                    const int& V_gcw,
                    const int& V_depth,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1
#if (NDIM == 3)
                    ,
                    const int& ilower2,
                    const int& iupper2
#endif
                    );

void C_TO_C_ANISO_F_LAPLACE_FC(double* F,
                               const int& F_gcw,
                               const double* alpha0,
                               const double* alpha1,
#if (NDIM == 3)
                               const double* alpha2,
#endif
                               const int& alpha_gcw,
                               const double* U,
                               const int& U_gcw,
                               const int& ilower0,
                               const int& iupper0,
                               const int& ilower1,
                               const int& iupper1,
#if (NDIM == 3)
                               const int& ilower2,
                               const int& iupper2,
#endif
                               const double* dx);

void C_TO_C_ANISO_F_LAPLACE_ADD_FC(double* F,
                                   const int& F_gcw,
                                   const double* alpha0,
                                   const double* alpha1,
#if (NDIM == 3)
                                   const double* alpha2,
#endif
                                   const int& alpha_gcw,
                                   const double* U,
                                   const int& U_gcw,
                                   const double& beta,
                                   const double* V,
                                   const int& V_gcw,
                                   const int& ilower0,
                                   const int& iupper0,
                                   const int& ilower1,
                                   const int& iupper1,
#if (NDIM == 3)
                                   const int& ilower2,
                                   const int& iupper2,
#endif
                                   const double* dx);

void C_TO_C_ANISO_F_DAMPED_LAPLACE_FC(double* F,
                                      const int& F_gcw,
                                      const double* alpha0,
                                      const double* alpha1,
#if (NDIM == 3)
                                      const double* alpha2,
#endif
                                      const int& alpha_gcw,
                                      const double& beta,
                                      const double* U,
                                      const int& U_gcw,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const double* dx);

void C_TO_C_ANISO_F_DAMPED_LAPLACE_ADD_FC(double* F,
                                          const int& F_gcw,
                                          const double* alpha0,
                                          const double* alpha1,
#if (NDIM == 3)
                                          const double* alpha2,
#endif
                                          const int& alpha_gcw,
                                          const double& beta,
                                          const double* U,
                                          const int& U_gcw,
                                          const double& gamma,
                                          const double* V,
                                          const int& V_gcw,
                                          const int& ilower0,
                                          const int& iupper0,
                                          const int& ilower1,
                                          const int& iupper1,
#if (NDIM == 3)
                                          const int& ilower2,
                                          const int& iupper2,
#endif
                                          const double* dx);

void C_TO_C_ANISO_S_LAPLACE_FC(double* F,
                               const int& F_gcw,
                               const double* alpha0,
                               const double* alpha1,
#if (NDIM == 3)
                               const double* alpha2,
#endif
                               const int& alpha_gcw,
                               const double* U,
                               const int& U_gcw,
                               const int& ilower0,
                               const int& iupper0,
                               const int& ilower1,
                               const int& iupper1,
#if (NDIM == 3)
                               const int& ilower2,
                               const int& iupper2,
#endif
                               const double* dx);

void C_TO_C_ANISO_S_LAPLACE_ADD_FC(double* F,
                                   const int& F_gcw,
                                   const double* alpha0,
                                   const double* alpha1,
#if (NDIM == 3)
                                   const double* alpha2,
#endif
                                   const int& alpha_gcw,
                                   const double* U,
                                   const int& U_gcw,
                                   const double& beta,
                                   const double* V,
                                   const int& V_gcw,
                                   const int& ilower0,
                                   const int& iupper0,
                                   const int& ilower1,
                                   const int& iupper1,
#if (NDIM == 3)
                                   const int& ilower2,
                                   const int& iupper2,
#endif
                                   const double* dx);

void C_TO_C_ANISO_S_DAMPED_LAPLACE_FC(double* F,
                                      const int& F_gcw,
                                      const double* alpha0,
                                      const double* alpha1,
#if (NDIM == 3)
                                      const double* alpha2,
#endif
                                      const int& alpha_gcw,
                                      const double& beta,
                                      const double* U,
                                      const int& U_gcw,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const double* dx);

void C_TO_C_ANISO_S_DAMPED_LAPLACE_ADD_FC(double* F,
                                          const int& F_gcw,
                                          const double* alpha0,
                                          const double* alpha1,
#if (NDIM == 3)
                                          const double* alpha2,
#endif
                                          const int& alpha_gcw,
                                          const double& beta,
                                          const double* U,
                                          const int& U_gcw,
                                          const double& gamma,
                                          const double* V,
                                          const int& V_gcw,
                                          const int& ilower0,
                                          const int& iupper0,
                                          const int& ilower1,
                                          const int& iupper1,
#if (NDIM == 3)
                                          const int& ilower2,
                                          const int& iupper2,
#endif
                                          const double* dx);

void C_TO_F_GRAD_FC(double* g0,
                    double* g1,
#if (NDIM == 3)
                    double* g2,
#endif
                    const int& g_gcw,
                    const double& alpha,
                    const double* U,
                    const int& U_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
#if (NDIM == 3)
                    const int& ilower2,
                    const int& iupper2,
#endif
                    const double* dx);

void C_TO_F_FLUX_FC(double* g0,
                    double* g1,
#if (NDIM == 3)
                    double* g2,
#endif
                    const int& g_gcw,
                    const double* alpha0,
                    const double* alpha1,
#if (NDIM == 3)
                    const double* alpha2,
#endif
                    const int& alpha_gcw,
                    const double* U,
                    const int& U_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
#if (NDIM == 3)
                    const int& ilower2,
                    const int& iupper2,
#endif
                    const double* dx);

void C_TO_F_ANISO_FLUX_FC(double* g0,
                          double* g1,
#if (NDIM == 3)
                          double* g2,
#endif
                          const int& g_gcw,
                          const double* alpha0,
                          const double* alpha1,
#if (NDIM == 3)
                          const double* alpha2,
#endif
                          const int& alpha_gcw,
                          const double* U,
                          const int& U_gcw,
                          const int& ilower0,
                          const int& iupper0,
                          const int& ilower1,
                          const int& iupper1,
#if (NDIM == 3)
                          const int& ilower2,
                          const int& iupper2,
#endif
                          const double* dx);

void C_TO_F_GRAD_ADD_FC(double* g0,
                        double* g1,
#if (NDIM == 3)
                        double* g2,
#endif
                        const int& g_gcw,
                        const double& alpha,
                        const double* U,
                        const int& U_gcw,
                        const double& beta,
                        const double* v0,
                        const double* v1,
#if (NDIM == 3)
                        const double* v2,
#endif
                        const int& v_gcw,
                        const int& ilower0,
                        const int& iupper0,
                        const int& ilower1,
                        const int& iupper1,
#if (NDIM == 3)
                        const int& ilower2,
                        const int& iupper2,
#endif
                        const double* dx);

void C_TO_F_INTERP_FC(double* u0,
                      double* u1,
#if (NDIM == 3)
                      double* u2,
#endif
                      const int& u_gcw,
                      const double* V,
                      const int& V_gcw,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1
#if (NDIM == 3)
                      ,
                      const int& ilower2,
                      const int& iupper2
#endif
                      );

void C_TO_S_GRAD_FC(double* g0,
                    double* g1,
#if (NDIM == 3)
                    double* g2,
#endif
                    const int& g_gcw,
                    const double& alpha,
                    const double* U,
                    const int& U_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
#if (NDIM == 3)
                    const int& ilower2,
                    const int& iupper2,
#endif
                    const double* dx);

void C_TO_S_FLUX_FC(double* g0,
                    double* g1,
#if (NDIM == 3)
                    double* g2,
#endif
                    const int& g_gcw,
                    const double* alpha0,
                    const double* alpha1,
#if (NDIM == 3)
                    const double* alpha2,
#endif
                    const int& alpha_gcw,
                    const double* U,
                    const int& U_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
#if (NDIM == 3)
                    const int& ilower2,
                    const int& iupper2,
#endif
                    const double* dx);

void C_TO_S_ANISO_FLUX_FC(double* g0,
                          double* g1,
#if (NDIM == 3)
                          double* g2,
#endif
                          const int& g_gcw,
                          const double* alpha0,
                          const double* alpha1,
#if (NDIM == 3)
                          const double* alpha2,
#endif
                          const int& alpha_gcw,
                          const double* U,
                          const int& U_gcw,
                          const int& ilower0,
                          const int& iupper0,
                          const int& ilower1,
                          const int& iupper1,
#if (NDIM == 3)
                          const int& ilower2,
                          const int& iupper2,
#endif
                          const double* dx);

void C_TO_S_GRAD_ADD_FC(double* g0,
                        double* g1,
#if (NDIM == 3)
                        double* g2,
#endif
                        const int& g_gcw,
                        const double& alpha,
                        const double* U,
                        const int& U_gcw,
                        const double& beta,
                        const double* v0,
                        const double* v1,
#if (NDIM == 3)
                        const double* v2,
#endif
                        const int& v_gcw,
                        const int& ilower0,
                        const int& iupper0,
                        const int& ilower1,
                        const int& iupper1,
#if (NDIM == 3)
                        const int& ilower2,
                        const int& iupper2,
#endif
                        const double* dx);

void C_TO_S_INTERP_FC(double* u0,
                      double* u1,
#if (NDIM == 3)
                      double* u2,
#endif
                      const int& u_gcw,
                      const double* V,
                      const int& V_gcw,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1
#if (NDIM == 3)
                      ,
                      const int& ilower2,
                      const int& iupper2
#endif
                      );

void F_TO_C_CURL_FC(double* W,
                    const int& W_gcw,
                    const double* u0,
                    const double* u1,
#if (NDIM == 3)
                    const double* u2,
#endif
                    const int& u_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
#if (NDIM == 3)
                    const int& ilower2,
                    const int& iupper2,
#endif
                    const double* dx);

void F_TO_C_DIV_FC(double* D,
                   const int& D_gcw,
                   const double& alpha,
                   const double* u0,
                   const double* u1,
#if (NDIM == 3)
                   const double* u2,
#endif
                   const int& u_gcw,
                   const int& ilower0,
                   const int& iupper0,
                   const int& ilower1,
                   const int& iupper1,
#if (NDIM == 3)
                   const int& ilower2,
                   const int& iupper2,
#endif
                   const double* dx);

void F_TO_C_DIV_ADD_FC(double* D,
                       const int& D_gcw,
                       const double& alpha,
                       const double* u0,
                       const double* u1,
#if (NDIM == 3)
                       const double* u2,
#endif
                       const int& u_gcw,
                       const double& beta,
                       const double* V,
                       const int& V_gcw,
                       const int& ilower0,
                       const int& iupper0,
                       const int& ilower1,
                       const int& iupper1,
#if (NDIM == 3)
                       const int& ilower2,
                       const int& iupper2,
#endif
                       const double* dx);

void F_TO_C_INTERP_FC(double* U,
                      const int& U_gcw,
                      const double* v0,
                      const double* v1,
#if (NDIM == 3)
                      const double* v2,
#endif
                      const int& v_gcw,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1
#if (NDIM == 3)
                      ,
                      const int& ilower2,
                      const int& iupper2
#endif
                      );

#if (NDIM == 3)
void F_TO_F_CURL_FC(double* w0,
                    double* w1,
                    double* w2,
                    const int& w_gcw,
                    const double* u0,
                    const double* u1,
                    const double* u2,
                    const int& u_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
                    const int& ilower2,
                    const int& iupper2,
                    const double* dx);
#endif

void S_TO_C_CURL_FC(double* W,
                    const int& W_gcw,
                    const double* u0,
                    const double* u1,
#if (NDIM == 3)
                    const double* u2,
#endif
                    const int& u_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
#if (NDIM == 3)
                    const int& ilower2,
                    const int& iupper2,
#endif
                    const double* dx);

void S_TO_C_DIV_FC(double* D,
                   const int& D_gcw,
                   const double& alpha,
                   const double* u0,
                   const double* u1,
#if (NDIM == 3)
                   const double* u2,
#endif
                   const int& u_gcw,
                   const int& ilower0,
                   const int& iupper0,
                   const int& ilower1,
                   const int& iupper1,
#if (NDIM == 3)
                   const int& ilower2,
                   const int& iupper2,
#endif
                   const double* dx);

void S_TO_C_DIV_ADD_FC(double* D,
                       const int& D_gcw,
                       const double& alpha,
                       const double* u0,
                       const double* u1,
#if (NDIM == 3)
                       const double* u2,
#endif
                       const int& u_gcw,
                       const double& beta,
                       const double* V,
                       const int& V_gcw,
                       const int& ilower0,
                       const int& iupper0,
                       const int& ilower1,
                       const int& iupper1,
#if (NDIM == 3)
                       const int& ilower2,
                       const int& iupper2,
#endif
                       const double* dx);

void S_TO_C_INTERP_FC(double* U,
                      const int& U_gcw,
                      const double* v0,
                      const double* v1,
#if (NDIM == 3)
                      const double* v2,
#endif
                      const int& v_gcw,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1
#if (NDIM == 3)
                      ,
                      const int& ilower2,
                      const int& iupper2
#endif
                      );

#if (NDIM == 3)
void S_TO_S_CURL_FC(double* w0,
                    double* w1,
                    double* w2,
                    const int& w_gcw,
                    const double* u0,
                    const double* u1,
                    const double* u2,
                    const int& u_gcw,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
                    const int& ilower2,
                    const int& iupper2,
                    const double* dx);
#endif
#if (NDIM == 2)
void S_TO_S_VC_LAPLACE_FC(double* f0,
                          double* f1,
#if (NDIM == 3)
                          double* f2,
#endif
                          const int& f_gcw,
                          const double& alpha,
                          const double& beta,
                          const double* mu,
                          const int& mu_gcw,
                          const double* u0,
                          const double* u1,
#if (NDIM == 3)
                          const double* u2,
#endif
                          const int& u_gcw,
                          const double& gamma,
                          const double* v0,
                          const double* v1,
#if (NDIM == 3)
                          const double* v2,
#endif
                          const int& v_gcw,
                          const int& ilower0,
                          const int& iupper0,
                          const int& ilower1,
                          const int& iupper1,
#if (NDIM == 3)
                          const int& ilower2,
                          const int& iupper2,
#endif
                          const double* dx);
#endif

#if (NDIM == 2)
void N_TO_S_ROT_FC(double* w0,
                   double* w1,
                   const int& w_ghosts,
                   const double* u0,
                   const int& u_ghosts,
                   const int& ilower0,
                   const int& iupper0,
                   const int& ilower1,
                   const int& iupper1,
                   const double* dx);

void C_TO_S_ROT_FC(double* w0,
                   double* w1,
                   const int& w_ghosts,
                   const double* u0,
                   const int& u_ghosts,
                   const int& ilower0,
                   const int& iupper0,
                   const int& ilower1,
                   const int& iupper1,
                   const double* dx);

void S_TO_N_CURL_FC(double* w0,
                    const int& w_ghosts,
                    const double* u0,
                    const double* u1,
                    const int& u_ghosts,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
                    const double* dx);
#endif

#if (NDIM == 3)
void E_TO_S_ROT_FC(double* w0,
                   double* w1,
                   double* w2,
                   const int& w_ghosts,
                   const double* u0,
                   const double* u1,
                   const double* u2,
                   const int& u_ghosts,
                   const int& ilower0,
                   const int& iupper0,
                   const int& ilower1,
                   const int& iupper1,
                   const int& ilower2,
                   const int& iupper2,
                   const double* dx);

void S_TO_E_CURL_FC(double* w0,
                    double* w1,
                    double* w2,
                    const int& w_ghosts,
                    const double* u0,
                    const double* u1,
                    const double* u2,
                    const int& u_ghosts,
                    const int& ilower0,
                    const int& iupper0,
                    const int& ilower1,
                    const int& iupper1,
                    const int& ilower2,
                    const int& iupper2,
                    const double* dx);
#endif
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PatchMathOps::PatchMathOps()
{
    // intentionally blank
    return;
} // PatchMathOps

PatchMathOps::~PatchMathOps()
{
    // intentionally blank
    return;
} // ~PatchMathOps

void PatchMathOps::curl(Pointer<CellData<NDIM, double> > dst,
                        const Pointer<CellData<NDIM, double> > src,
                        const Pointer<Patch<NDIM> > patch) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const W = dst->getPointer();
    const int W_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src->getPointer();
    const int U_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (W_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (src == dst)
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src == dst." << std::endl);
    }

    const Box<NDIM>& U_box = src->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src has insufficient ghost cell width" << std::endl);
    }

    const int W_depth = dst->getDepth();

    if (
#if (NDIM == 2)
        (W_depth != 1)
#endif
#if (NDIM == 3)
            (W_depth != NDIM)
#endif
                )
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst has incorrect depth" << std::endl);
    }

    const int U_depth = src->getDepth();

    if (U_depth != NDIM)
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src has incorrect depth" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    C_TO_C_CURL_FC(W, W_ghosts, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                   patch_box.upper(1),
#if (NDIM == 3)
                   patch_box.lower(2), patch_box.upper(2),
#endif
                   dx);
    return;
} // curl

void PatchMathOps::curl(Pointer<CellData<NDIM, double> > dst,
                        const Pointer<FaceData<NDIM, double> > src,
                        const Pointer<Patch<NDIM> > patch) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const W = dst->getPointer();
    const int W_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src->getPointer(0);
    const double* const u1 = src->getPointer(1);
#if (NDIM == 3)
    const double* const u2 = src->getPointer(2);
#endif
    const int u_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (W_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (src == dst)
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src == dst." << std::endl);
    }

    const Box<NDIM>& U_box = src->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src has insufficient ghost cell width" << std::endl);
    }

    const int W_depth = dst->getDepth();

    if (
#if (NDIM == 2)
        (W_depth != 1)
#endif
#if (NDIM == 3)
            (W_depth != NDIM)
#endif
                )
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst has incorrect depth" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    F_TO_C_CURL_FC(W, W_ghosts, u0, u1,
#if (NDIM == 3)
                   u2,
#endif
                   u_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                   patch_box.lower(2), patch_box.upper(2),
#endif
                   dx);
    return;
} // curl

void PatchMathOps::curl(Pointer<FaceData<NDIM, double> > dst,
                        const Pointer<FaceData<NDIM, double> > src,
                        const Pointer<Patch<NDIM> > patch) const
{
#if (NDIM != 3)
    TBOX_ERROR("PatchMathOps::curl():\n"
               << "  not implemented for NDIM != 3" << std::endl);
    NULL_USE(dst);
    NULL_USE(src);
    NULL_USE(patch);
#endif
#if (NDIM == 3)
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const w0 = dst->getPointer(0);
    double* const w1 = dst->getPointer(1);
    double* const w2 = dst->getPointer(2);
    const int w_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src->getPointer(0);
    const double* const u1 = src->getPointer(1);
    const double* const u2 = src->getPointer(2);
    const int u_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (w_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (src == dst)
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src == dst." << std::endl);
    }

    const Box<NDIM>& U_box = src->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    F_TO_F_CURL_FC(w0, w1, w2, w_ghosts, u0, u1, u2, u_ghosts, patch_box.lower(0), patch_box.upper(0),
                   patch_box.lower(1), patch_box.upper(1), patch_box.lower(2), patch_box.upper(2), dx);
#endif
    return;
} // curl

void PatchMathOps::curl(Pointer<CellData<NDIM, double> > dst,
                        const Pointer<SideData<NDIM, double> > src,
                        const Pointer<Patch<NDIM> > patch) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const W = dst->getPointer();
    const int W_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src->getPointer(0);
    const double* const u1 = src->getPointer(1);
#if (NDIM == 3)
    const double* const u2 = src->getPointer(2);
#endif
    const int u_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (W_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (src == dst)
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src == dst." << std::endl);
    }

    const Box<NDIM>& U_box = src->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src has insufficient ghost cell width" << std::endl);
    }

    const int W_depth = dst->getDepth();

    if (
#if (NDIM == 2)
        (W_depth != 1)
#endif
#if (NDIM == 3)
            (W_depth != NDIM)
#endif
                )
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst has incorrect depth" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    S_TO_C_CURL_FC(W, W_ghosts, u0, u1,
#if (NDIM == 3)
                   u2,
#endif
                   u_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                   patch_box.lower(2), patch_box.upper(2),
#endif
                   dx);
    return;
} // curl

void PatchMathOps::curl(Pointer<SideData<NDIM, double> > dst,
                        const Pointer<SideData<NDIM, double> > src,
                        const Pointer<Patch<NDIM> > patch) const
{
#if (NDIM != 3)
    TBOX_ERROR("PatchMathOps::curl():\n"
               << "  not implemented for NDIM != 3" << std::endl);
    NULL_USE(dst);
    NULL_USE(src);
    NULL_USE(patch);
#endif
#if (NDIM == 3)
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const w0 = dst->getPointer(0);
    double* const w1 = dst->getPointer(1);
    double* const w2 = dst->getPointer(2);
    const int w_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src->getPointer(0);
    const double* const u1 = src->getPointer(1);
    const double* const u2 = src->getPointer(2);
    const int u_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (w_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (src == dst)
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src == dst." << std::endl);
    }

    const Box<NDIM>& U_box = src->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    S_TO_S_CURL_FC(w0, w1, w2, w_ghosts, u0, u1, u2, u_ghosts, patch_box.lower(0), patch_box.upper(0),
                   patch_box.lower(1), patch_box.upper(1), patch_box.lower(2), patch_box.upper(2), dx);
#endif
    return;
} // curl

void PatchMathOps::curl(Pointer<NodeData<NDIM, double> > dst,
                        const Pointer<SideData<NDIM, double> > src,
                        const Pointer<Patch<NDIM> > patch) const
{
#if (NDIM != 2)
    TBOX_ERROR("PatchMathOps::curl():\n"
               << "  not implemented for NDIM != 2" << std::endl);
    NULL_USE(dst);
    NULL_USE(src);
    NULL_USE(patch);
#endif
#if (NDIM == 2)
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const w0 = dst->getPointer(0);
    const int w_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src->getPointer(0);
    const double* const u1 = src->getPointer(1);
    const int u_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (w_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (src == dst)
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src == dst." << std::endl);
    }

    const Box<NDIM>& U_box = src->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    S_TO_N_CURL_FC(w0, w_ghosts, u0, u1, u_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                   patch_box.upper(1), dx);
#endif
    return;
} // curl

void PatchMathOps::curl(Pointer<EdgeData<NDIM, double> > dst,
                        const Pointer<SideData<NDIM, double> > src,
                        const Pointer<Patch<NDIM> > patch) const
{
#if (NDIM != 3)
    TBOX_ERROR("PatchMathOps::curl():\n"
               << "  not implemented for NDIM != 3" << std::endl);
    NULL_USE(dst);
    NULL_USE(src);
    NULL_USE(patch);
#endif
#if (NDIM == 3)
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const w0 = dst->getPointer(0);
    double* const w1 = dst->getPointer(1);
    double* const w2 = dst->getPointer(2);
    const int w_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src->getPointer(0);
    const double* const u1 = src->getPointer(1);
    const double* const u2 = src->getPointer(2);
    const int u_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (w_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (src == dst)
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src == dst." << std::endl);
    }

    const Box<NDIM>& U_box = src->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  src has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::curl():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    S_TO_E_CURL_FC(w0, w1, w2, w_ghosts, u0, u1, u2, u_ghosts, patch_box.lower(0), patch_box.upper(0),
                   patch_box.lower(1), patch_box.upper(1), patch_box.lower(2), patch_box.upper(2), dx);
#endif
    return;
} // curl

void PatchMathOps::rot(Pointer<SideData<NDIM, double> > dst,
                       const Pointer<NodeData<NDIM, double> > src,
                       const Pointer<Patch<NDIM> > patch,
                       CartSideRobinPhysBdryOp* bc_op,
                       const double fill_time) const
{
#if (NDIM != 2)
    TBOX_ERROR("PatchMathOps::rot():\n"
               << "  not implemented for NDIM != 2" << std::endl);
    NULL_USE(dst);
    NULL_USE(src);
    NULL_USE(patch);
    NULL_USE(bc_op);
    NULL_USE(fill_time);
#endif
#if (NDIM == 2)
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const w0 = dst->getPointer(0);
    double* const w1 = dst->getPointer(1);
    const int w_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src->getPointer(0);
    const int u_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (w_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (src == dst)
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  src == dst." << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    if (!bc_op)
    {
        N_TO_S_ROT_FC(w0, w1, w_ghosts, u0, u_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                      patch_box.upper(1), dx);
    }
    else
    {
        IntVector<NDIM> op_gcw = IntVector<NDIM>(1);

        IntVector<NDIM> u_gcw = IntVector<NDIM>::max(IntVector<NDIM>(2), src->getGhostCellWidth());
        NodeData<NDIM, double> u_data(patch_box, src->getDepth(), u_gcw);
        const double* const u0 = u_data.getPointer(0);
        const int u_ghosts = (u_data.getGhostCellWidth() - op_gcw).max();
        u_data.fillAll(0.0);
        Box<NDIM> copy_box = patch_box;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const int lower = 0;
            const int upper = 1;
            if (!pgeom->getTouchesRegularBoundary(axis, lower)) copy_box.lower(axis) -= u_gcw(axis);
            if (!pgeom->getTouchesRegularBoundary(axis, upper)) copy_box.upper(axis) += u_gcw(axis);
        }
        u_data.copyOnBox(*src, copy_box);

        IntVector<NDIM> w_gcw = IntVector<NDIM>::max(IntVector<NDIM>(2), dst->getGhostCellWidth());
        SideData<NDIM, double> w_data(patch_box, dst->getDepth(), w_gcw);
        double* const w0 = w_data.getPointer(0);
        double* const w1 = w_data.getPointer(1);
        const int w_ghosts = (w_data.getGhostCellWidth() - op_gcw).max();

        const Box<NDIM> op_box = Box<NDIM>::grow(patch_box, op_gcw);

        N_TO_S_ROT_FC(w0, w1, w_ghosts, u0, u_ghosts, op_box.lower(0), op_box.upper(0), op_box.lower(1),
                      op_box.upper(1), dx);

        dst->copyOnBox(w_data, op_box);

        bc_op->accumulateFromPhysicalBoundaryData(*patch, fill_time, op_gcw);
    }
#endif
    return;
} // rot

void PatchMathOps::rot(Pointer<SideData<NDIM, double> > dst,
                       const Pointer<CellData<NDIM, double> > src,
                       const Pointer<Patch<NDIM> > patch,
                       CartSideRobinPhysBdryOp* bc_op,
                       const double fill_time) const
{
#if (NDIM != 2)
    TBOX_ERROR("PatchMathOps::rot():\n"
               << "  not implemented for NDIM != 2" << std::endl);
    NULL_USE(dst);
    NULL_USE(src);
    NULL_USE(patch);
    NULL_USE(bc_op);
    NULL_USE(fill_time);
#endif
#if (NDIM == 2)
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const w0 = dst->getPointer(0);
    double* const w1 = dst->getPointer(1);
    const int w_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src->getPointer(0);
    const int u_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (w_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (src == dst)
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  src == dst." << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    const Box<NDIM>& U_box = src->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  src has insufficient ghost cell width" << std::endl);
    }
#endif

    if (!bc_op)
    {
        C_TO_S_ROT_FC(w0, w1, w_ghosts, u0, u_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                      patch_box.upper(1), dx);
    }
    else
    {
        IntVector<NDIM> op_gcw = IntVector<NDIM>(1);

        IntVector<NDIM> u_gcw = IntVector<NDIM>::max(IntVector<NDIM>(2), src->getGhostCellWidth());
        CellData<NDIM, double> u_data(patch_box, src->getDepth(), u_gcw);
        const double* const u0 = u_data.getPointer(0);
        const int u_ghosts = (u_data.getGhostCellWidth() - op_gcw).max();
        u_data.fillAll(0.0);
        Box<NDIM> copy_box = patch_box;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const int lower = 0;
            const int upper = 1;
            if (!pgeom->getTouchesRegularBoundary(axis, lower)) copy_box.lower(axis) -= u_gcw(axis);
            if (!pgeom->getTouchesRegularBoundary(axis, upper)) copy_box.upper(axis) += u_gcw(axis);
        }
        u_data.copyOnBox(*src, copy_box);

        IntVector<NDIM> w_gcw = IntVector<NDIM>::max(IntVector<NDIM>(2), dst->getGhostCellWidth());
        SideData<NDIM, double> w_data(patch_box, dst->getDepth(), w_gcw);
        double* const w0 = w_data.getPointer(0);
        double* const w1 = w_data.getPointer(1);
        const int w_ghosts = (w_data.getGhostCellWidth() - op_gcw).max();

        const Box<NDIM> op_box = Box<NDIM>::grow(patch_box, op_gcw);

        C_TO_S_ROT_FC(w0, w1, w_ghosts, u0, u_ghosts, op_box.lower(0), op_box.upper(0), op_box.lower(1),
                      op_box.upper(1), dx);

        dst->copyOnBox(w_data, op_box);

        bc_op->accumulateFromPhysicalBoundaryData(*patch, fill_time, op_gcw);
    }
#endif
    return;
} // rot

void PatchMathOps::rot(Pointer<SideData<NDIM, double> > dst,
                       const Pointer<EdgeData<NDIM, double> > src,
                       const Pointer<Patch<NDIM> > patch,
                       CartSideRobinPhysBdryOp* bc_op,
                       const double fill_time) const
{
#if (NDIM != 3)
    TBOX_ERROR("PatchMathOps::rot():\n"
               << "  not implemented for NDIM != 3" << std::endl);
    NULL_USE(dst);
    NULL_USE(src);
    NULL_USE(patch);
    NULL_USE(bc_op);
    NULL_USE(fill_time);
#endif
#if (NDIM == 3)
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const w0 = dst->getPointer(0);
    double* const w1 = dst->getPointer(1);
    double* const w2 = dst->getPointer(2);
    const int w_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src->getPointer(0);
    const double* const u1 = src->getPointer(1);
    const double* const u2 = src->getPointer(2);
    const int u_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (w_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (src == dst)
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  src == dst." << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    if (!bc_op)
    {
        E_TO_S_ROT_FC(w0, w1, w2, w_ghosts, u0, u1, u2, u_ghosts, patch_box.lower(0), patch_box.upper(0),
                      patch_box.lower(1), patch_box.upper(1), patch_box.lower(2), patch_box.upper(2), dx);
    }
    else
    {
        IntVector<NDIM> op_gcw = IntVector<NDIM>(1);

        IntVector<NDIM> u_gcw = IntVector<NDIM>::max(IntVector<NDIM>(2), src->getGhostCellWidth());
        EdgeData<NDIM, double> u_data(patch_box, src->getDepth(), u_gcw);
        const double* const u0 = u_data.getPointer(0);
        const double* const u1 = u_data.getPointer(1);
        const double* const u2 = u_data.getPointer(2);
        const int u_ghosts = (u_data.getGhostCellWidth() - op_gcw).max();
        u_data.fillAll(0.0);
        Box<NDIM> copy_box = patch_box;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const int lower = 0;
            const int upper = 1;
            if (!pgeom->getTouchesRegularBoundary(axis, lower)) copy_box.lower(axis) -= u_gcw(axis);
            if (!pgeom->getTouchesRegularBoundary(axis, upper)) copy_box.upper(axis) += u_gcw(axis);
        }
        u_data.copyOnBox(*src, copy_box);

        IntVector<NDIM> w_gcw = IntVector<NDIM>::max(IntVector<NDIM>(2), dst->getGhostCellWidth());
        SideData<NDIM, double> w_data(patch_box, dst->getDepth(), w_gcw);
        double* const w0 = w_data.getPointer(0);
        double* const w1 = w_data.getPointer(1);
        double* const w2 = w_data.getPointer(1);
        const int w_ghosts = (w_data.getGhostCellWidth() - op_gcw).max();

        const Box<NDIM> op_box = Box<NDIM>::grow(patch_box, op_gcw);

        E_TO_S_ROT_FC(w0, w1, w2, w_ghosts, u0, u1, u2, u_ghosts, op_box.lower(0), op_box.upper(0), op_box.lower(1),
                      op_box.upper(1), op_box.lower(2), op_box.upper(2), dx);

        dst->copyOnBox(w_data, op_box);

        bc_op->accumulateFromPhysicalBoundaryData(*patch, fill_time, op_gcw);
    }
#endif
    return;
} // rot

void PatchMathOps::rot(Pointer<SideData<NDIM, double> > dst,
                       const Pointer<SideData<NDIM, double> > src,
                       const Pointer<Patch<NDIM> > patch,
                       CartSideRobinPhysBdryOp* bc_op,
                       const double fill_time) const
{
    if (!bc_op)
    {
        this->curl(dst, src, patch);
    }
    else
    {
#if (NDIM != 3)
        TBOX_ERROR("PatchMathOps::rot():\n"
                   << "  not implemented for NDIM != 3" << std::endl);
        NULL_USE(dst);
        NULL_USE(src);
        NULL_USE(patch);
        NULL_USE(bc_op);
        NULL_USE(fill_time);
#endif
#if (NDIM == 3)
        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        const Box<NDIM>& patch_box = patch->getBox();

        IntVector<NDIM> op_gcw = IntVector<NDIM>(1);

        IntVector<NDIM> u_gcw = IntVector<NDIM>::max(IntVector<NDIM>(2), src->getGhostCellWidth());
        SideData<NDIM, double> u_data(patch_box, src->getDepth(), u_gcw);
        const double* const u0 = src->getPointer(0);
        const double* const u1 = src->getPointer(1);
        const double* const u2 = src->getPointer(2);
        const int u_ghosts = (src->getGhostCellWidth()).max();

        IntVector<NDIM> w_gcw = IntVector<NDIM>::max(IntVector<NDIM>(2), dst->getGhostCellWidth());
        SideData<NDIM, double> w_data(patch_box, dst->getDepth(), w_gcw);
        double* const w0 = w_data.getPointer(0);
        double* const w1 = w_data.getPointer(1);
        double* const w2 = w_data.getPointer(1);
        const int w_ghosts = (w_data.getGhostCellWidth() - op_gcw).max();

        const Box<NDIM> op_box = Box<NDIM>::grow(patch_box, op_gcw);

#if !defined(NDEBUG)
        if (w_ghosts != (dst->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::rot():\n"
                       << "  dst does not have uniform ghost cell widths" << std::endl);
        }

        if (u_ghosts != (src->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::rot():\n"
                       << "  src does not have uniform ghost cell widths" << std::endl);
        }

        if (src == dst)
        {
            TBOX_ERROR("PatchMathOps::rot():\n"
                       << "  src == dst." << std::endl);
        }

        const Box<NDIM>& U_box = src->getGhostBox();
        const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

        if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
        {
            TBOX_ERROR("PatchMathOps::rot():\n"
                       << "  src has insufficient ghost cell width" << std::endl);
        }

        if (patch_box != dst->getBox())
        {
            TBOX_ERROR("PatchMathOps::rot():\n"
                       << "  dst and src must live on the same patch" << std::endl);
        }

        if (patch_box != src->getBox())
        {
            TBOX_ERROR("PatchMathOps::rot():\n"
                       << "  dst and src must live on the same patch" << std::endl);
        }
#endif

        S_TO_S_CURL_FC(w0, w1, w2, w_ghosts, u0, u1, u2, u_ghosts, op_box.lower(0), op_box.upper(0), op_box.lower(1),
                       op_box.upper(1), op_box.lower(2), op_box.upper(2), dx);

        dst->copyOnBox(w_data, op_box);

        bc_op->accumulateFromPhysicalBoundaryData(*patch, fill_time, op_gcw);
#endif
    }
    return;
} // rot

void PatchMathOps::div(Pointer<CellData<NDIM, double> > dst,
                       const double alpha,
                       const Pointer<CellData<NDIM, double> > src1,
                       const double beta,
                       const Pointer<CellData<NDIM, double> > src2,
                       const Pointer<Patch<NDIM> > patch,
                       const int l,
                       const int m) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const D = dst->getPointer(l);
    const int D_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer();
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (D_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (src1 == dst)
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  src1 == dst." << std::endl);
    }

    if ((src1 == src2) && (beta != 0.0))
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  src1 == src2 but beta is nonzero." << std::endl);
    }

    const Box<NDIM>& U_box = src1->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    const int U_depth = src1->getDepth();

    if (U_depth != NDIM)
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  src1 has incorrect depth" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (beta == 0.0))
    {
        C_TO_C_DIV_FC(D, D_ghosts, alpha, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                      patch_box.upper(1),
#if (NDIM == 3)
                      patch_box.lower(2), patch_box.upper(2),
#endif
                      dx);
    }
    else
    {
        const double* const V = src2->getPointer(m);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::div():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::div():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        C_TO_C_DIV_ADD_FC(D, D_ghosts, alpha, U, U_ghosts, beta, V, V_ghosts, patch_box.lower(0), patch_box.upper(0),
                          patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                          patch_box.lower(2), patch_box.upper(2),
#endif
                          dx);
    }
    return;
} // div

void PatchMathOps::div(Pointer<CellData<NDIM, double> > dst,
                       const double alpha,
                       const Pointer<FaceData<NDIM, double> > src1,
                       const double beta,
                       const Pointer<CellData<NDIM, double> > src2,
                       const Pointer<Patch<NDIM> > patch,
                       const int l,
                       const int m) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const D = dst->getPointer(l);
    const int D_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src1->getPointer(0);
    const double* const u1 = src1->getPointer(1);
#if (NDIM == 3)
    const double* const u2 = src1->getPointer(2);
#endif
    const int u_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (D_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (beta == 0.0))
    {
        F_TO_C_DIV_FC(D, D_ghosts, alpha, u0, u1,
#if (NDIM == 3)
                      u2,
#endif
                      u_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                      patch_box.lower(2), patch_box.upper(2),
#endif
                      dx);
    }
    else
    {
        const double* const V = src2->getPointer(m);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::div():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::div():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        F_TO_C_DIV_ADD_FC(D, D_ghosts, alpha, u0, u1,
#if (NDIM == 3)
                          u2,
#endif
                          u_ghosts, beta, V, V_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                          patch_box.upper(1),
#if (NDIM == 3)
                          patch_box.lower(2), patch_box.upper(2),
#endif
                          dx);
    }
    return;
} // div

void PatchMathOps::div(Pointer<CellData<NDIM, double> > dst,
                       const double alpha,
                       const Pointer<SideData<NDIM, double> > src1,
                       const double beta,
                       const Pointer<CellData<NDIM, double> > src2,
                       const Pointer<Patch<NDIM> > patch,
                       const int l,
                       const int m) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const D = dst->getPointer(l);
    const int D_ghosts = (dst->getGhostCellWidth()).max();

    const double* const u0 = src1->getPointer(0);
    const double* const u1 = src1->getPointer(1);
#if (NDIM == 3)
    const double* const u2 = src1->getPointer(2);
#endif
    const int u_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (D_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::div():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (beta == 0.0))
    {
        S_TO_C_DIV_FC(D, D_ghosts, alpha, u0, u1,
#if (NDIM == 3)
                      u2,
#endif
                      u_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                      patch_box.lower(2), patch_box.upper(2),
#endif
                      dx);
    }
    else
    {
        const double* const V = src2->getPointer(m);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::div():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::div():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        S_TO_C_DIV_ADD_FC(D, D_ghosts, alpha, u0, u1,
#if (NDIM == 3)
                          u2,
#endif
                          u_ghosts, beta, V, V_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                          patch_box.upper(1),
#if (NDIM == 3)
                          patch_box.lower(2), patch_box.upper(2),
#endif
                          dx);
    }
    return;
} // div

void PatchMathOps::grad(Pointer<CellData<NDIM, double> > dst,
                        const double alpha,
                        const Pointer<CellData<NDIM, double> > src1,
                        const double beta,
                        const Pointer<CellData<NDIM, double> > src2,
                        const Pointer<Patch<NDIM> > patch,
                        const int l) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const G = dst->getPointer();
    const int G_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(l);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (G_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    const int G_depth = dst->getDepth();

    if (G_depth != NDIM)
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst has incorrect depth" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (src1 == dst)
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 == dst." << std::endl);
    }

    if ((src1 == src2) && (beta != 0.0))
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 == src2 but beta is nonzero." << std::endl);
    }

    const Box<NDIM>& U_box = src1->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (beta == 0.0))
    {
        C_TO_C_GRAD_FC(G, G_ghosts, alpha, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                       patch_box.upper(1),
#if (NDIM == 3)
                       patch_box.lower(2), patch_box.upper(2),
#endif
                       dx);
    }
    else
    {
        const double* const V = src2->getPointer();
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::grad():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        const int V_depth = src2->getDepth();

        if (V_depth != NDIM)
        {
            TBOX_ERROR("PatchMathOps::grad():\n"
                       << "  src2 has incorrect depth" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::grad():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        C_TO_C_GRAD_ADD_FC(G, G_ghosts, alpha, U, U_ghosts, beta, V, V_ghosts, patch_box.lower(0), patch_box.upper(0),
                           patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                           patch_box.lower(2), patch_box.upper(2),
#endif
                           dx);
    }
    return;
} // grad

void PatchMathOps::grad(Pointer<FaceData<NDIM, double> > dst,
                        const double alpha,
                        const Pointer<CellData<NDIM, double> > src1,
                        const double beta,
                        const Pointer<FaceData<NDIM, double> > src2,
                        const Pointer<Patch<NDIM> > patch,
                        const int l) const
{
    // Compute the gradient.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const g0 = dst->getPointer(0);
    double* const g1 = dst->getPointer(1);
#if (NDIM == 3)
    double* const g2 = dst->getPointer(2);
#endif
    const int g_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(l);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (g_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    const Box<NDIM>& U_box = src1->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (beta == 0.0))
    {
        C_TO_F_GRAD_FC(g0, g1,
#if (NDIM == 3)
                       g2,
#endif
                       g_ghosts, alpha, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                       patch_box.upper(1),
#if (NDIM == 3)
                       patch_box.lower(2), patch_box.upper(2),
#endif
                       dx);
    }
    else
    {
        const double* const v0 = src2->getPointer(0);
        const double* const v1 = src2->getPointer(1);
#if (NDIM == 3)
        const double* const v2 = src2->getPointer(2);
#endif
        const int v_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (v_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::grad():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::grad():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        C_TO_F_GRAD_ADD_FC(g0, g1,
#if (NDIM == 3)
                           g2,
#endif
                           g_ghosts, alpha, U, U_ghosts, beta, v0, v1,
#if (NDIM == 3)
                           v2,
#endif
                           v_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                           patch_box.lower(2), patch_box.upper(2),
#endif
                           dx);
    }
    return;
} // grad

void PatchMathOps::grad(Pointer<SideData<NDIM, double> > dst,
                        const double alpha,
                        const Pointer<CellData<NDIM, double> > src1,
                        const double beta,
                        const Pointer<SideData<NDIM, double> > src2,
                        const Pointer<Patch<NDIM> > patch,
                        const int l) const
{
    // Compute the gradient.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const g0 = dst->getPointer(0);
    double* const g1 = dst->getPointer(1);
#if (NDIM == 3)
    double* const g2 = dst->getPointer(2);
#endif
    const int g_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(l);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (g_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    const Box<NDIM>& U_box = src1->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (beta == 0.0))
    {
        C_TO_S_GRAD_FC(g0, g1,
#if (NDIM == 3)
                       g2,
#endif
                       g_ghosts, alpha, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                       patch_box.upper(1),
#if (NDIM == 3)
                       patch_box.lower(2), patch_box.upper(2),
#endif
                       dx);
    }
    else
    {
        const double* const v0 = src2->getPointer(0);
        const double* const v1 = src2->getPointer(1);
#if (NDIM == 3)
        const double* const v2 = src2->getPointer(2);
#endif
        const int v_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (v_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::grad():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::grad():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        C_TO_S_GRAD_ADD_FC(g0, g1,
#if (NDIM == 3)
                           g2,
#endif
                           g_ghosts, alpha, U, U_ghosts, beta, v0, v1,
#if (NDIM == 3)
                           v2,
#endif
                           v_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                           patch_box.lower(2), patch_box.upper(2),
#endif
                           dx);
    }
    return;
} // grad

void PatchMathOps::grad(Pointer<FaceData<NDIM, double> > dst,
                        const Pointer<FaceData<NDIM, double> > alpha,
                        const Pointer<CellData<NDIM, double> > src1,
                        const double beta,
                        const Pointer<FaceData<NDIM, double> > src2,
                        const Pointer<Patch<NDIM> > patch,
                        const int l) const
{
    // Compute the gradient.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const g0 = dst->getPointer(0);
    double* const g1 = dst->getPointer(1);
#if (NDIM == 3)
    double* const g2 = dst->getPointer(2);
#endif
    const int g_ghosts = (dst->getGhostCellWidth()).max();

    const double* const alpha0 = alpha->getPointer(0);
    const double* const alpha1 = alpha->getPointer(1);
#if (NDIM == 3)
    const double* const alpha2 = alpha->getPointer(2);
#endif
    const int alpha_ghosts = (alpha->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(l);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (g_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    const int g_depth = dst->getDepth();

    if (g_depth != 1)
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst has incorrect depth" << std::endl);
    }

    if (alpha_ghosts != (alpha->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  alpha does not have uniform ghost cell widths" << std::endl);
    }

    const int a_depth = alpha->getDepth();

    if (a_depth != 1 && a_depth != NDIM)
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  alpha has incorrect depth" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    const Box<NDIM>& U_box = src1->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != alpha->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif
    if (alpha->getDepth() == 1)
    {
        C_TO_F_FLUX_FC(g0, g1,
#if (NDIM == 3)
                       g2,
#endif
                       g_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                       alpha2,
#endif
                       alpha_ghosts, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                       patch_box.upper(1),
#if (NDIM == 3)
                       patch_box.lower(2), patch_box.upper(2),
#endif
                       dx);

        if (src2 && (beta != 0.0))
        {
#if !defined(NDEBUG)
            const int v_ghosts = (src2->getGhostCellWidth()).max();

            if (v_ghosts != (src2->getGhostCellWidth()).min())
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  src2 does not have uniform ghost cell widths" << std::endl);
            }

            const int v_depth = dst->getDepth();

            if (v_depth != 1)
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  src2 has incorrect depth" << std::endl);
            }

            if (patch_box != src2->getBox())
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
            }
#endif
            PatchFaceDataOpsReal<NDIM, double> patch_fc_data_ops;
            patch_fc_data_ops.axpy(dst, beta, src2, dst, patch_box);
        }
    }
    else
    {
        C_TO_F_ANISO_FLUX_FC(g0, g1,
#if (NDIM == 3)
                             g2,
#endif
                             g_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                             alpha2,
#endif
                             alpha_ghosts, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                             patch_box.upper(1),
#if (NDIM == 3)
                             patch_box.lower(2), patch_box.upper(2),
#endif
                             dx);

        // Account for non-zero beta.
        if (src2 && (beta != 0.0))
        {
#if !defined(NDEBUG)
            const int v_ghosts = (src2->getGhostCellWidth()).max();

            if (v_ghosts != (src2->getGhostCellWidth()).min())
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  src2 does not have uniform ghost cell widths" << std::endl);
            }

            const int v_depth = dst->getDepth();

            if (v_depth != 1)
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  src2 has incorrect depth" << std::endl);
            }

            if (patch_box != src2->getBox())
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
            }
#endif
            PatchFaceDataOpsReal<NDIM, double> patch_fc_data_ops;
            patch_fc_data_ops.axpy(dst, beta, src2, dst, patch_box);
        }
    }
    return;
} // grad

void PatchMathOps::grad(Pointer<SideData<NDIM, double> > dst,
                        const Pointer<SideData<NDIM, double> > alpha,
                        const Pointer<CellData<NDIM, double> > src1,
                        const double beta,
                        const Pointer<SideData<NDIM, double> > src2,
                        const Pointer<Patch<NDIM> > patch,
                        const int l) const
{
    // Compute the gradient.
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const g0 = dst->getPointer(0);
    double* const g1 = dst->getPointer(1);
#if (NDIM == 3)
    double* const g2 = dst->getPointer(2);
#endif
    const int g_ghosts = (dst->getGhostCellWidth()).max();

    const double* const alpha0 = alpha->getPointer(0);
    const double* const alpha1 = alpha->getPointer(1);
#if (NDIM == 3)
    const double* const alpha2 = alpha->getPointer(2);
#endif
    const int alpha_ghosts = (alpha->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(l);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (g_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    const int g_depth = dst->getDepth();

    if (g_depth != 1)
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst has incorrect depth" << std::endl);
    }

    if (alpha_ghosts != (alpha->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  alpha does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    const Box<NDIM>& U_box = src1->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != alpha->getBox())
    {
        TBOX_ERROR("PatchMathOps::grad():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif
    if (alpha->getDepth() == 1)
    {
        C_TO_S_FLUX_FC(g0, g1,
#if (NDIM == 3)
                       g2,
#endif
                       g_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                       alpha2,
#endif
                       alpha_ghosts, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                       patch_box.upper(1),
#if (NDIM == 3)
                       patch_box.lower(2), patch_box.upper(2),
#endif
                       dx);

        if (src2 && (beta != 0.0))
        {
#if !defined(NDEBUG)
            const int v_ghosts = (src2->getGhostCellWidth()).max();

            if (v_ghosts != (src2->getGhostCellWidth()).min())
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  src2 does not have uniform ghost cell widths" << std::endl);
            }

            const int v_depth = dst->getDepth();

            if (v_depth != 1)
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  src2 has incorrect depth" << std::endl);
            }

            if (patch_box != src2->getBox())
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
            }
#endif
            PatchSideDataOpsReal<NDIM, double> patch_sc_data_ops;
            patch_sc_data_ops.axpy(dst, beta, src2, dst, patch_box);
        }
    }
    else
    {
        C_TO_S_ANISO_FLUX_FC(g0, g1,
#if (NDIM == 3)
                             g2,
#endif
                             g_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                             alpha2,
#endif
                             alpha_ghosts, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                             patch_box.upper(1),
#if (NDIM == 3)
                             patch_box.lower(2), patch_box.upper(2),
#endif
                             dx);

        // Account for non-zero beta.
        if (src2 && (beta != 0.0))
        {
#if !defined(NDEBUG)
            const int v_ghosts = (src2->getGhostCellWidth()).max();

            if (v_ghosts != (src2->getGhostCellWidth()).min())
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  src2 does not have uniform ghost cell widths" << std::endl);
            }

            const int v_depth = dst->getDepth();

            if (v_depth != 1)
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  src2 has incorrect depth" << std::endl);
            }

            if (patch_box != src2->getBox())
            {
                TBOX_ERROR("PatchMathOps::grad():\n"
                           << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
            }
#endif
            PatchSideDataOpsReal<NDIM, double> patch_sc_data_ops;
            patch_sc_data_ops.axpy(dst, beta, src2, dst, patch_box);
        }
    }
    return;
} // grad

void PatchMathOps::interp(Pointer<CellData<NDIM, double> > dst,
                          const Pointer<FaceData<NDIM, double> > src,
                          const Pointer<Patch<NDIM> > patch) const
{
    const int U_ghosts = (dst->getGhostCellWidth()).max();
    const int v_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (dst->getDepth() != NDIM * src->getDepth())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  src and dst have incompatible depths" << std::endl);
    }

    if (U_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (v_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    for (int depth = 0; depth < src->getDepth(); ++depth)
    {
        // Interpolate.
        double* const U = dst->getPointer(NDIM * depth);

        const double* const v0 = src->getPointer(0, depth);
        const double* const v1 = src->getPointer(1, depth);
#if (NDIM == 3)
        const double* const v2 = src->getPointer(2, depth);
#endif

        F_TO_C_INTERP_FC(U, U_ghosts, v0, v1,
#if (NDIM == 3)
                         v2,
#endif
                         v_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1), patch_box.upper(1)
#if (NDIM == 3)
                                                                                                   ,
                         patch_box.lower(2), patch_box.upper(2)
#endif
                             );
    }
    return;
} // interp

void PatchMathOps::interp(Pointer<CellData<NDIM, double> > dst,
                          const Pointer<SideData<NDIM, double> > src,
                          const Pointer<Patch<NDIM> > patch) const
{
    const int U_ghosts = (dst->getGhostCellWidth()).max();
    const int v_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (dst->getDepth() != NDIM * src->getDepth())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  src and dst have incompatible depths" << std::endl);
    }

    if (U_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (v_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    for (int depth = 0; depth < src->getDepth(); ++depth)
    {
        // Interpolate.
        double* const U = dst->getPointer(NDIM * depth);

        const double* const v0 = src->getPointer(0, depth);
        const double* const v1 = src->getPointer(1, depth);
#if (NDIM == 3)
        const double* const v2 = src->getPointer(2, depth);
#endif

        S_TO_C_INTERP_FC(U, U_ghosts, v0, v1,
#if (NDIM == 3)
                         v2,
#endif
                         v_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1), patch_box.upper(1)
#if (NDIM == 3)
                                                                                                   ,
                         patch_box.lower(2), patch_box.upper(2)
#endif
                             );
    }
    return;
} // interp

void PatchMathOps::interp(Pointer<FaceData<NDIM, double> > dst,
                          const Pointer<CellData<NDIM, double> > src,
                          const Pointer<Patch<NDIM> > patch) const
{
    const int u_ghosts = (dst->getGhostCellWidth()).max();
    const int V_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (NDIM * dst->getDepth() != src->getDepth())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  src and dst have incompatible depths" << std::endl);
    }

    if (u_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (V_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    const Box<NDIM>& V_box = src->getGhostBox();
    const Box<NDIM> V_box_shrunk = Box<NDIM>::grow(V_box, -1);

    if ((!V_box_shrunk.contains(patch_box.lower())) || (!V_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  src has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    for (int depth = 0; depth < dst->getDepth(); ++depth)
    {
        // Interpolate.
        double* const u0 = dst->getPointer(0, depth);
        double* const u1 = dst->getPointer(1, depth);
#if (NDIM == 3)
        double* const u2 = dst->getPointer(2, depth);
#endif
        const double* const V = src->getPointer(NDIM * depth);

        C_TO_F_INTERP_FC(u0, u1,
#if (NDIM == 3)
                         u2,
#endif
                         u_ghosts, V, V_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                         patch_box.upper(1)
#if (NDIM == 3)
                             ,
                         patch_box.lower(2), patch_box.upper(2)
#endif
                             );
    }
    return;
} // interp

void PatchMathOps::interp(Pointer<SideData<NDIM, double> > dst,
                          const Pointer<CellData<NDIM, double> > src,
                          const Pointer<Patch<NDIM> > patch) const
{
    const int u_ghosts = (dst->getGhostCellWidth()).max();
    const int V_ghosts = (src->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (NDIM * dst->getDepth() != src->getDepth())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  src and dst have incompatible depths" << std::endl);
    }

    if (u_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (V_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    const Box<NDIM>& V_box = src->getGhostBox();
    const Box<NDIM> V_box_shrunk = Box<NDIM>::grow(V_box, -1);

    if ((!V_box_shrunk.contains(patch_box.lower())) || (!V_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  src has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::interp():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    for (int depth = 0; depth < dst->getDepth(); ++depth)
    {
        // Interpolate.
        double* const u0 = dst->getPointer(0, depth);
        double* const u1 = dst->getPointer(1, depth);
#if (NDIM == 3)
        double* const u2 = dst->getPointer(2, depth);
#endif
        const double* const V = src->getPointer(NDIM * depth);

        C_TO_S_INTERP_FC(u0, u1,
#if (NDIM == 3)
                         u2,
#endif
                         u_ghosts, V, V_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                         patch_box.upper(1)
#if (NDIM == 3)
                             ,
                         patch_box.lower(2), patch_box.upper(2)
#endif
                             );
    }
    return;
} // interp

void PatchMathOps::laplace(Pointer<CellData<NDIM, double> > dst,
                           const double alpha,
                           const double beta,
                           const Pointer<CellData<NDIM, double> > src1,
                           const double gamma,
                           const Pointer<CellData<NDIM, double> > src2,
                           const Pointer<Patch<NDIM> > patch,
                           const int l,
                           const int m,
                           const int n) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const F = dst->getPointer(l);
    const int F_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(m);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (F_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (src1 == dst)
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 == dst." << std::endl);
    }

    if ((src1 == src2) && (gamma != 0.0))
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 == src2 but gamma is nonzero." << std::endl);
    }

    const Box<NDIM>& U_box = src1->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (gamma == 0.0))
    {
        if (beta == 0.0)
        {
            LAPLACE_FC(F, F_ghosts, alpha, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                       patch_box.upper(1),
#if (NDIM == 3)
                       patch_box.lower(2), patch_box.upper(2),
#endif
                       dx);
        }
        else
        {
            DAMPED_LAPLACE_FC(F, F_ghosts, alpha, beta, U, U_ghosts, patch_box.lower(0), patch_box.upper(0),
                              patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                              patch_box.lower(2), patch_box.upper(2),
#endif
                              dx);
        }
    }
    else
    {
        const double* const V = src2->getPointer(n);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::laplace():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::laplace():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        if (beta == 0.0)
        {
            LAPLACE_ADD_FC(F, F_ghosts, alpha, U, U_ghosts, gamma, V, V_ghosts, patch_box.lower(0), patch_box.upper(0),
                           patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                           patch_box.lower(2), patch_box.upper(2),
#endif
                           dx);
        }
        else
        {
            DAMPED_LAPLACE_ADD_FC(F, F_ghosts, alpha, beta, U, U_ghosts, gamma, V, V_ghosts, patch_box.lower(0),
                                  patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                                  patch_box.lower(2), patch_box.upper(2),
#endif
                                  dx);
        }
    }
    return;
} // laplace

void PatchMathOps::laplace(Pointer<SideData<NDIM, double> > dst,
                           const double alpha,
                           const double beta,
                           const Pointer<SideData<NDIM, double> > src1,
                           const double gamma,
                           const Pointer<SideData<NDIM, double> > src2,
                           const Pointer<Patch<NDIM> > patch,
                           const int l,
                           const int m,
                           const int n) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    boost::array<double*, NDIM> F;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F[d] = dst->getPointer(d, l);
    }
    const int F_ghosts = (dst->getGhostCellWidth()).max();

    boost::array<const double*, NDIM> U;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        U[d] = src1->getPointer(d, m);
    }
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (F_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (src1 == dst)
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 == dst." << std::endl);
    }

    if ((src1 == src2) && (gamma != 0.0))
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 == src2 but gamma is nonzero." << std::endl);
    }

    const Box<NDIM>& U_box = src1->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (gamma == 0.0))
    {
        boost::array<int, NDIM> ilower, iupper;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            for (unsigned int dd = 0; dd < NDIM; ++dd)
            {
                ilower[dd] = patch_box.lower(dd);
                iupper[dd] = patch_box.upper(dd);
            }
            iupper[d] += 1;

            if (beta == 0.0)
            {
                LAPLACE_FC(F[d], F_ghosts, alpha, U[d], U_ghosts, ilower[0], iupper[0], ilower[1], iupper[1],
#if (NDIM == 3)
                           ilower[2], iupper[2],
#endif
                           dx);
            }
            else
            {
                DAMPED_LAPLACE_FC(F[d], F_ghosts, alpha, beta, U[d], U_ghosts, ilower[0], iupper[0], ilower[1],
                                  iupper[1],
#if (NDIM == 3)
                                  ilower[2], iupper[2],
#endif
                                  dx);
            }
        }
    }
    else
    {
        boost::array<const double*, NDIM> V;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            V[d] = src2->getPointer(d, n);
        }
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::laplace():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::laplace():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        boost::array<int, NDIM> ilower, iupper;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            for (unsigned int dd = 0; dd < NDIM; ++dd)
            {
                ilower[dd] = patch_box.lower(dd);
                iupper[dd] = patch_box.upper(dd);
            }
            iupper[d] += 1;

            if (beta == 0.0)
            {
                LAPLACE_ADD_FC(F[d], F_ghosts, alpha, U[d], U_ghosts, gamma, V[d], V_ghosts, ilower[0], iupper[0],
                               ilower[1], iupper[1],
#if (NDIM == 3)
                               ilower[2], iupper[2],
#endif
                               dx);
            }
            else
            {
                DAMPED_LAPLACE_ADD_FC(F[d], F_ghosts, alpha, beta, U[d], U_ghosts, gamma, V[d], V_ghosts, ilower[0],
                                      iupper[0], ilower[1], iupper[1],
#if (NDIM == 3)
                                      ilower[2], iupper[2],
#endif
                                      dx);
            }
        }
    }
    return;
} // laplace

void PatchMathOps::laplace(Pointer<CellData<NDIM, double> > dst,
                           const Pointer<FaceData<NDIM, double> > alpha,
                           const double beta,
                           const Pointer<CellData<NDIM, double> > src1,
                           const double gamma,
                           const Pointer<CellData<NDIM, double> > src2,
                           const Pointer<Patch<NDIM> > patch,
                           const int l,
                           const int m,
                           const int n) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const F = dst->getPointer(l);
    const int F_ghosts = (dst->getGhostCellWidth()).max();

    const double* const alpha0 = alpha->getPointer(0);
    const double* const alpha1 = alpha->getPointer(1);
#if (NDIM == 3)
    const double* const alpha2 = alpha->getPointer(2);
#endif
    const int alpha_ghosts = (alpha->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(m);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (F_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (alpha_ghosts != (alpha->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  alpha does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (src1 == dst)
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 == dst." << std::endl);
    }

    if ((src1 == src2) && (gamma != 0.0))
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 == src2 but gamma is nonzero." << std::endl);
    }

    const Box<NDIM>& U_box = src1->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != alpha->getBox())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (alpha->getDepth() != 1)
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  alpha has invalid depth" << std::endl);
    }
#endif

    if (!src2 || (gamma == 0.0))
    {
        if (beta == 0.0)
        {
            C_TO_C_ANISO_F_LAPLACE_FC(F, F_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                                      alpha2,
#endif
                                      alpha_ghosts, U, U_ghosts, patch_box.lower(0), patch_box.upper(0),
                                      patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                                      patch_box.lower(2), patch_box.upper(2),
#endif
                                      dx);
        }
        else
        {
            C_TO_C_ANISO_F_DAMPED_LAPLACE_FC(F, F_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                                             alpha2,
#endif
                                             alpha_ghosts, beta, U, U_ghosts, patch_box.lower(0), patch_box.upper(0),
                                             patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                                             patch_box.lower(2), patch_box.upper(2),
#endif
                                             dx);
        }
    }
    else
    {
        const double* const V = src2->getPointer(n);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::laplace():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::laplace():\n"
                       << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        if (beta == 0.0)
        {
            C_TO_C_ANISO_F_LAPLACE_ADD_FC(F, F_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                                          alpha2,
#endif
                                          alpha_ghosts, U, U_ghosts, gamma, V, V_ghosts, patch_box.lower(0),
                                          patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                                          patch_box.lower(2), patch_box.upper(2),
#endif
                                          dx);
        }
        else
        {
            C_TO_C_ANISO_F_DAMPED_LAPLACE_ADD_FC(F, F_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                                                 alpha2,
#endif
                                                 alpha_ghosts, beta, U, U_ghosts, gamma, V, V_ghosts,
                                                 patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                                                 patch_box.upper(1),
#if (NDIM == 3)
                                                 patch_box.lower(2), patch_box.upper(2),
#endif
                                                 dx);
        }
    }
    return;
} // laplace

void PatchMathOps::laplace(Pointer<CellData<NDIM, double> > dst,
                           const Pointer<SideData<NDIM, double> > alpha,
                           const double beta,
                           const Pointer<CellData<NDIM, double> > src1,
                           const double gamma,
                           const Pointer<CellData<NDIM, double> > src2,
                           const Pointer<Patch<NDIM> > patch,
                           const int l,
                           const int m,
                           const int n) const
{
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const F = dst->getPointer(l);
    const int F_ghosts = (dst->getGhostCellWidth()).max();

    const double* const alpha0 = alpha->getPointer(0);
    const double* const alpha1 = alpha->getPointer(1);
#if (NDIM == 3)
    const double* const alpha2 = alpha->getPointer(2);
#endif
    const int alpha_ghosts = (alpha->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(m);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (F_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (alpha_ghosts != (alpha->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  alpha does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (src1 == dst)
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 == dst." << std::endl);
    }

    if ((src1 == src2) && (gamma != 0.0))
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 == src2 but gamma is nonzero." << std::endl);
    }

    const Box<NDIM>& U_box = src1->getGhostBox();
    const Box<NDIM> U_box_shrunk = Box<NDIM>::grow(U_box, -1);

    if ((!U_box_shrunk.contains(patch_box.lower())) || (!U_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != alpha->getBox())
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (alpha->getDepth() != 1)
    {
        TBOX_ERROR("PatchMathOps::laplace():\n"
                   << "  alpha has invalid depth" << std::endl);
    }
#endif

    if (!src2 || (gamma == 0.0))
    {
        if (beta == 0.0)
        {
            C_TO_C_ANISO_S_LAPLACE_FC(F, F_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                                      alpha2,
#endif
                                      alpha_ghosts, U, U_ghosts, patch_box.lower(0), patch_box.upper(0),
                                      patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                                      patch_box.lower(2), patch_box.upper(2),
#endif
                                      dx);
        }
        else
        {
            C_TO_C_ANISO_S_DAMPED_LAPLACE_FC(F, F_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                                             alpha2,
#endif
                                             alpha_ghosts, beta, U, U_ghosts, patch_box.lower(0), patch_box.upper(0),
                                             patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                                             patch_box.lower(2), patch_box.upper(2),
#endif
                                             dx);
        }
    }
    else
    {
        const double* const V = src2->getPointer(n);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::laplace():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::laplace():\n"
                       << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        if (beta == 0.0)
        {
            C_TO_C_ANISO_S_LAPLACE_ADD_FC(F, F_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                                          alpha2,
#endif
                                          alpha_ghosts, U, U_ghosts, gamma, V, V_ghosts, patch_box.lower(0),
                                          patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                                          patch_box.lower(2), patch_box.upper(2),
#endif
                                          dx);
        }
        else
        {
            C_TO_C_ANISO_S_DAMPED_LAPLACE_ADD_FC(F, F_ghosts, alpha0, alpha1,
#if (NDIM == 3)
                                                 alpha2,
#endif
                                                 alpha_ghosts, beta, U, U_ghosts, gamma, V, V_ghosts,
                                                 patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                                                 patch_box.upper(1),
#if (NDIM == 3)
                                                 patch_box.lower(2), patch_box.upper(2),
#endif
                                                 dx);
        }
    }
    return;
} // laplace

void PatchMathOps::vc_laplace(Pointer<SideData<NDIM, double> > dst,
                              const double alpha,
                              const double beta,
                              const Pointer<NodeData<NDIM, double> > coef,
                              const Pointer<SideData<NDIM, double> > src1,
                              const double gamma_in,
                              const Pointer<SideData<NDIM, double> > src2_in,
                              const Pointer<Patch<NDIM> > patch,
                              const int l,
                              const int m,
                              const int n) const
{
#if (NDIM == 2)
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    double* const f0 = dst->getPointer(0, l);
    double* const f1 = dst->getPointer(1, l);
#if (NDIM == 3)
    double* const f2 = dst->getPointer(2, l);
#endif
    const int f_ghosts = (dst->getGhostCellWidth()).max();

    const double* const mu = coef->getPointer();
    const int mu_ghosts = (coef->getGhostCellWidth()).max();

    const double* const u0 = src1->getPointer(0, m);
    const double* const u1 = src1->getPointer(1, m);
#if (NDIM == 3)
    const double* const u2 = src1->getPointer(2, m);
#endif
    const int u_ghosts = (src1->getGhostCellWidth()).max();

    const double gamma = (src2_in ? gamma_in : 0.0);
    const Pointer<SideData<NDIM, double> > src2 = (src2_in ? src2_in : src1);
    const double* const v0 = src2->getPointer(0, n);
    const double* const v1 = src2->getPointer(1, n);
#if (NDIM == 3)
    const double* const v2 = src2->getPointer(2, n);
#endif
    const int v_ghosts = (src2->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (f_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (mu_ghosts != (coef->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  coef does not have uniform ghost cell widths" << std::endl);
    }

    if (u_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (v_ghosts != (src2->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  src2 does not have uniform ghost cell widths" << std::endl);
    }

    if (src1 == dst)
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  src1 == dst." << std::endl);
    }

    if ((src1 == src2) && (gamma != 0.0))
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  src1 == src2 but gamma is nonzero." << std::endl);
    }

    const Box<NDIM>& mu_box = coef->getGhostBox();
    const Box<NDIM> mu_box_shrunk = Box<NDIM>::grow(mu_box, -1);

    if ((!mu_box_shrunk.contains(patch_box.lower())) || (!mu_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  coef has insufficient ghost cell width" << std::endl);
    }

    const Box<NDIM>& u_box = src1->getGhostBox();
    const Box<NDIM> u_box_shrunk = Box<NDIM>::grow(u_box, -1);

    if ((!u_box_shrunk.contains(patch_box.lower())) || (!u_box_shrunk.contains(patch_box.upper())))
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  src1 has insufficient ghost cell width" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  dst, coef, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != coef->getBox())
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  dst, coef, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  dst, coef, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src2->getBox())
    {
        TBOX_ERROR("PatchMathOps::vc_laplace():\n"
                   << "  dst, coef, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    S_TO_S_VC_LAPLACE_FC(f0, f1,
#if (NDIM == 3)
                         f2,
#endif
                         f_ghosts, alpha, beta, mu, mu_ghosts, u0, u1,
#if (NDIM == 3)
                         u2,
#endif
                         u_ghosts, gamma, v0, v1,
#if (NDIM == 3)
                         v2,
#endif
                         v_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                         patch_box.lower(2), patch_box.upper(2),
#endif
                         dx);
#endif
#if (NDIM == 3)
    TBOX_ERROR("PatchMathOps::vc_laplace():\n"
               << "  not presently implemented for NDIM = 3." << std::endl);
    NULL_USE(dst);
    NULL_USE(alpha);
    NULL_USE(beta);
    NULL_USE(coef);
    NULL_USE(src1);
    NULL_USE(gamma_in);
    NULL_USE(src2_in);
    NULL_USE(patch);
    NULL_USE(l);
    NULL_USE(m);
    NULL_USE(n);
#endif
    return;
} // vc_laplace

void PatchMathOps::pointwiseMultiply(Pointer<CellData<NDIM, double> > dst,
                                     const double alpha,
                                     const Pointer<CellData<NDIM, double> > src1,
                                     const double beta,
                                     const Pointer<CellData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k) const
{
    double* const D = dst->getPointer(i);
    const int D_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(j);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (D_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (beta == 0.0))
    {
        MULTIPLY1_FC(D, D_ghosts, alpha, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                     patch_box.upper(1)
#if (NDIM == 3)
                         ,
                     patch_box.lower(2), patch_box.upper(2)
#endif
                         );
    }
    else
    {
        const double* const V = src2->getPointer(k);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        MULTIPLY_ADD1_FC(D, D_ghosts, alpha, U, U_ghosts, beta, V, V_ghosts, patch_box.lower(0), patch_box.upper(0),
                         patch_box.lower(1), patch_box.upper(1)
#if (NDIM == 3)
                                                 ,
                         patch_box.lower(2), patch_box.upper(2)
#endif
                             );
    }
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<CellData<NDIM, double> > dst,
                                     const Pointer<CellData<NDIM, double> > alpha,
                                     const Pointer<CellData<NDIM, double> > src1,
                                     const double beta,
                                     const Pointer<CellData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k,
                                     const int l) const
{
    double* const D = dst->getPointer(i);
    const int D_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(j);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const double* const A = alpha->getPointer(l);
    const int A_ghosts = (alpha->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (D_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (A_ghosts != (alpha->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  alpha does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != alpha->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (beta == 0.0))
    {
        MULTIPLY2_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                     patch_box.upper(1)
#if (NDIM == 3)
                         ,
                     patch_box.lower(2), patch_box.upper(2)
#endif
                         );
    }
    else
    {
        const double* const V = src2->getPointer(k);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        MULTIPLY_ADD2_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, beta, V, V_ghosts, patch_box.lower(0),
                         patch_box.upper(0), patch_box.lower(1), patch_box.upper(1)
#if (NDIM == 3)
                                                                     ,
                         patch_box.lower(2), patch_box.upper(2)
#endif
                             );
    }
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<CellData<NDIM, double> > dst,
                                     const Pointer<CellData<NDIM, double> > alpha,
                                     const Pointer<CellData<NDIM, double> > src1,
                                     const Pointer<CellData<NDIM, double> > beta,
                                     const Pointer<CellData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k,
                                     const int l,
                                     const int m) const
{
    if (!src2)
    {
        const double zero = 0.0;
        pointwiseMultiply(dst, alpha, src1, zero, src2, patch, i, j, k, l);
        return;
    }

    double* const D = dst->getPointer(i);
    const int D_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(j);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const double* const V = src2->getPointer(k);
    const int V_ghosts = (src2->getGhostCellWidth()).max();

    const double* const A = alpha->getPointer(l);
    const int A_ghosts = (alpha->getGhostCellWidth()).max();

    const double* const B = beta->getPointer(m);
    const int B_ghosts = (beta->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (D_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (V_ghosts != (src2->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  src2 does not have uniform ghost cell widths" << std::endl);
    }

    if (A_ghosts != (alpha->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  alpha does not have uniform ghost cell widths" << std::endl);
    }

    if (B_ghosts != (beta->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  beta does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src2->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != alpha->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != beta->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
    }
#endif

    MULTIPLY_ADD3_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, B, B_ghosts, V, V_ghosts, patch_box.lower(0),
                     patch_box.upper(0), patch_box.lower(1), patch_box.upper(1)
#if (NDIM == 3)
                                                                 ,
                     patch_box.lower(2), patch_box.upper(2)
#endif
                         );
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<FaceData<NDIM, double> > dst,
                                     const double alpha,
                                     const Pointer<FaceData<NDIM, double> > src1,
                                     const double beta,
                                     const Pointer<FaceData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k) const
{
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        double* const D = dst->getPointer(axis, i);
        const int D_ghosts = (dst->getGhostCellWidth()).max();

        const double* const U = src1->getPointer(axis, j);
        const int U_ghosts = (src1->getGhostCellWidth()).max();

        const Box<NDIM>& patch_box = patch->getBox();
        const Box<NDIM> data_box = FaceGeometry<NDIM>::toFaceBox(patch_box, axis);

#if !defined(NDEBUG)
        if (D_ghosts != (dst->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst does not have uniform ghost cell widths" << std::endl);
        }

        if (U_ghosts != (src1->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src1 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != dst->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != src1->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif

        if (!src2 || (beta == 0.0))
        {
            MULTIPLY1_FC(D, D_ghosts, alpha, U, U_ghosts, data_box.lower(0), data_box.upper(0), data_box.lower(1),
                         data_box.upper(1)
#if (NDIM == 3)
                             ,
                         data_box.lower(2), data_box.upper(2)
#endif
                             );
        }
        else
        {
            const double* const V = src2->getPointer(axis, k);
            const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
            if (V_ghosts != (src2->getGhostCellWidth()).min())
            {
                TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                           << "  src2 does not have uniform ghost cell widths" << std::endl);
            }

            if (patch_box != src2->getBox())
            {
                TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                           << "  dst, src1, and src2 must all live on the same patch" << std::endl);
            }
#endif
            MULTIPLY_ADD1_FC(D, D_ghosts, alpha, U, U_ghosts, beta, V, V_ghosts, data_box.lower(0), data_box.upper(0),
                             data_box.lower(1), data_box.upper(1)
#if (NDIM == 3)
                                                    ,
                             data_box.lower(2), data_box.upper(2)
#endif
                                 );
        }
    }
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<FaceData<NDIM, double> > dst,
                                     const Pointer<FaceData<NDIM, double> > alpha,
                                     const Pointer<FaceData<NDIM, double> > src1,
                                     const double beta,
                                     const Pointer<FaceData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k,
                                     const int l) const
{
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        double* const D = dst->getPointer(axis, i);
        const int D_ghosts = (dst->getGhostCellWidth()).max();

        const double* const U = src1->getPointer(axis, j);
        const int U_ghosts = (src1->getGhostCellWidth()).max();

        const double* const A = alpha->getPointer(axis, l);
        const int A_ghosts = (alpha->getGhostCellWidth()).max();

        const Box<NDIM>& patch_box = patch->getBox();
        const Box<NDIM> data_box = FaceGeometry<NDIM>::toFaceBox(patch_box, axis);

#if !defined(NDEBUG)
        if (D_ghosts != (dst->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst does not have uniform ghost cell widths" << std::endl);
        }

        if (U_ghosts != (src1->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src1 does not have uniform ghost cell widths" << std::endl);
        }

        if (A_ghosts != (alpha->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  alpha does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != dst->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != src1->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != alpha->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif

        if (!src2 || (beta == 0.0))
        {
            MULTIPLY2_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, data_box.lower(0), data_box.upper(0), data_box.lower(1),
                         data_box.upper(1)
#if (NDIM == 3)
                             ,
                         data_box.lower(2), data_box.upper(2)
#endif
                             );
        }
        else
        {
            const double* const V = src2->getPointer(axis, k);
            const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
            if (V_ghosts != (src2->getGhostCellWidth()).min())
            {
                TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                           << "  src2 does not have uniform ghost cell widths" << std::endl);
            }

            if (patch_box != src2->getBox())
            {
                TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                           << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
            }
#endif
            MULTIPLY_ADD2_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, beta, V, V_ghosts, data_box.lower(0),
                             data_box.upper(0), data_box.lower(1), data_box.upper(1)
#if (NDIM == 3)
                                                                       ,
                             data_box.lower(2), data_box.upper(2)
#endif
                                 );
        }
    }
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<FaceData<NDIM, double> > dst,
                                     const Pointer<FaceData<NDIM, double> > alpha,
                                     const Pointer<FaceData<NDIM, double> > src1,
                                     const Pointer<FaceData<NDIM, double> > beta,
                                     const Pointer<FaceData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k,
                                     const int l,
                                     const int m) const
{
    if (!src2)
    {
        const double zero = 0.0;
        pointwiseMultiply(dst, alpha, src1, zero, src2, patch, i, j, k, l);
        return;
    }

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        double* const D = dst->getPointer(axis, i);
        const int D_ghosts = (dst->getGhostCellWidth()).max();

        const double* const U = src1->getPointer(axis, j);
        const int U_ghosts = (src1->getGhostCellWidth()).max();

        const double* const V = src2->getPointer(axis, k);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

        const double* const A = alpha->getPointer(axis, l);
        const int A_ghosts = (alpha->getGhostCellWidth()).max();

        const double* const B = beta->getPointer(axis, m);
        const int B_ghosts = (beta->getGhostCellWidth()).max();

        const Box<NDIM>& patch_box = patch->getBox();
        const Box<NDIM> data_box = FaceGeometry<NDIM>::toFaceBox(patch_box, axis);

#if !defined(NDEBUG)
        if (D_ghosts != (dst->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst does not have uniform ghost cell widths" << std::endl);
        }

        if (U_ghosts != (src1->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src1 does not have uniform ghost cell widths" << std::endl);
        }

        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (A_ghosts != (alpha->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  alpha does not have uniform ghost cell widths" << std::endl);
        }

        if (B_ghosts != (beta->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  beta does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != dst->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != src1->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != alpha->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != beta->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
        }
#endif

        MULTIPLY_ADD3_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, B, B_ghosts, V, V_ghosts, data_box.lower(0),
                         data_box.upper(0), data_box.lower(1), data_box.upper(1)
#if (NDIM == 3)
                                                                   ,
                         data_box.lower(2), data_box.upper(2)
#endif
                             );
    }
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<NodeData<NDIM, double> > dst,
                                     const double alpha,
                                     const Pointer<NodeData<NDIM, double> > src1,
                                     const double beta,
                                     const Pointer<NodeData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k) const
{
    double* const D = dst->getPointer(i);
    const int D_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(j);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();
    const Box<NDIM> data_box = NodeGeometry<NDIM>::toNodeBox(patch_box);

#if !defined(NDEBUG)
    if (D_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (beta == 0.0))
    {
        MULTIPLY1_FC(D, D_ghosts, alpha, U, U_ghosts, data_box.lower(0), data_box.upper(0), data_box.lower(1),
                     data_box.upper(1)
#if (NDIM == 3)
                         ,
                     data_box.lower(2), data_box.upper(2)
#endif
                         );
    }
    else
    {
        const double* const V = src2->getPointer(k);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        MULTIPLY_ADD1_FC(D, D_ghosts, alpha, U, U_ghosts, beta, V, V_ghosts, data_box.lower(0), data_box.upper(0),
                         data_box.lower(1), data_box.upper(1)
#if (NDIM == 3)
                                                ,
                         data_box.lower(2), data_box.upper(2)
#endif
                             );
    }
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<NodeData<NDIM, double> > dst,
                                     const Pointer<NodeData<NDIM, double> > alpha,
                                     const Pointer<NodeData<NDIM, double> > src1,
                                     const double beta,
                                     const Pointer<NodeData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k,
                                     const int l) const
{
    double* const D = dst->getPointer(i);
    const int D_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(j);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const double* const A = alpha->getPointer(l);
    const int A_ghosts = (alpha->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();
    const Box<NDIM> data_box = NodeGeometry<NDIM>::toNodeBox(patch_box);

#if !defined(NDEBUG)
    if (D_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (A_ghosts != (alpha->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  alpha does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != alpha->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
    }
#endif

    if (!src2 || (beta == 0.0))
    {
        MULTIPLY2_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, data_box.lower(0), data_box.upper(0), data_box.lower(1),
                     data_box.upper(1)
#if (NDIM == 3)
                         ,
                     data_box.lower(2), data_box.upper(2)
#endif
                         );
    }
    else
    {
        const double* const V = src2->getPointer(k);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif
        MULTIPLY_ADD2_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, beta, V, V_ghosts, data_box.lower(0), data_box.upper(0),
                         data_box.lower(1), data_box.upper(1)
#if (NDIM == 3)
                                                ,
                         data_box.lower(2), data_box.upper(2)
#endif
                             );
    }
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<NodeData<NDIM, double> > dst,
                                     const Pointer<NodeData<NDIM, double> > alpha,
                                     const Pointer<NodeData<NDIM, double> > src1,
                                     const Pointer<NodeData<NDIM, double> > beta,
                                     const Pointer<NodeData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k,
                                     const int l,
                                     const int m) const
{
    if (!src2)
    {
        const double zero = 0.0;
        pointwiseMultiply(dst, alpha, src1, zero, src2, patch, i, j, k, l);
        return;
    }

    double* const D = dst->getPointer(i);
    const int D_ghosts = (dst->getGhostCellWidth()).max();

    const double* const U = src1->getPointer(j);
    const int U_ghosts = (src1->getGhostCellWidth()).max();

    const double* const V = src2->getPointer(k);
    const int V_ghosts = (src2->getGhostCellWidth()).max();

    const double* const A = alpha->getPointer(l);
    const int A_ghosts = (alpha->getGhostCellWidth()).max();

    const double* const B = beta->getPointer(m);
    const int B_ghosts = (beta->getGhostCellWidth()).max();

    const Box<NDIM>& patch_box = patch->getBox();
    const Box<NDIM> data_box = NodeGeometry<NDIM>::toNodeBox(patch_box);

#if !defined(NDEBUG)
    if (D_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    if (U_ghosts != (src1->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  src1 does not have uniform ghost cell widths" << std::endl);
    }

    if (V_ghosts != (src2->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  src2 does not have uniform ghost cell widths" << std::endl);
    }

    if (A_ghosts != (alpha->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  alpha does not have uniform ghost cell widths" << std::endl);
    }

    if (B_ghosts != (beta->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  beta does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src1->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != src2->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != alpha->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
    }

    if (patch_box != beta->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                   << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
    }
#endif

    MULTIPLY_ADD3_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, B, B_ghosts, V, V_ghosts, data_box.lower(0),
                     data_box.upper(0), data_box.lower(1), data_box.upper(1)
#if (NDIM == 3)
                                                               ,
                     data_box.lower(2), data_box.upper(2)
#endif
                         );
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<SideData<NDIM, double> > dst,
                                     const double alpha,
                                     const Pointer<SideData<NDIM, double> > src1,
                                     const double beta,
                                     const Pointer<SideData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k) const
{
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        double* const D = dst->getPointer(axis, i);
        const int D_ghosts = (dst->getGhostCellWidth()).max();

        const double* const U = src1->getPointer(axis, j);
        const int U_ghosts = (src1->getGhostCellWidth()).max();

        const Box<NDIM>& patch_box = patch->getBox();
        const Box<NDIM> data_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);

#if !defined(NDEBUG)
        if (D_ghosts != (dst->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst does not have uniform ghost cell widths" << std::endl);
        }

        if (U_ghosts != (src1->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src1 does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != dst->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != src1->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif

        if (!src2 || (beta == 0.0))
        {
            MULTIPLY1_FC(D, D_ghosts, alpha, U, U_ghosts, data_box.lower(0), data_box.upper(0), data_box.lower(1),
                         data_box.upper(1)
#if (NDIM == 3)
                             ,
                         data_box.lower(2), data_box.upper(2)
#endif
                             );
        }
        else
        {
            const double* const V = src2->getPointer(axis, k);
            const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
            if (V_ghosts != (src2->getGhostCellWidth()).min())
            {
                TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                           << "  src2 does not have uniform ghost cell widths" << std::endl);
            }

            if (patch_box != src2->getBox())
            {
                TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                           << "  dst, src1, and src2 must all live on the same patch" << std::endl);
            }
#endif
            MULTIPLY_ADD1_FC(D, D_ghosts, alpha, U, U_ghosts, beta, V, V_ghosts, data_box.lower(0), data_box.upper(0),
                             data_box.lower(1), data_box.upper(1)
#if (NDIM == 3)
                                                    ,
                             data_box.lower(2), data_box.upper(2)
#endif
                                 );
        }
    }
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<SideData<NDIM, double> > dst,
                                     const Pointer<SideData<NDIM, double> > alpha,
                                     const Pointer<SideData<NDIM, double> > src1,
                                     const double beta,
                                     const Pointer<SideData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k,
                                     const int l) const
{
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        double* const D = dst->getPointer(axis, i);
        const int D_ghosts = (dst->getGhostCellWidth()).max();

        const double* const U = src1->getPointer(axis, j);
        const int U_ghosts = (src1->getGhostCellWidth()).max();

        const double* const A = alpha->getPointer(axis, l);
        const int A_ghosts = (alpha->getGhostCellWidth()).max();

        const Box<NDIM>& patch_box = patch->getBox();
        const Box<NDIM> data_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);

#if !defined(NDEBUG)
        if (D_ghosts != (dst->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst does not have uniform ghost cell widths" << std::endl);
        }

        if (U_ghosts != (src1->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src1 does not have uniform ghost cell widths" << std::endl);
        }

        if (A_ghosts != (alpha->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  alpha does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != dst->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != src1->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != alpha->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
        }
#endif

        if (!src2 || (beta == 0.0))
        {
            MULTIPLY2_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, data_box.lower(0), data_box.upper(0), data_box.lower(1),
                         data_box.upper(1)
#if (NDIM == 3)
                             ,
                         data_box.lower(2), data_box.upper(2)
#endif
                             );
        }
        else
        {
            const double* const V = src2->getPointer(axis, k);
            const int V_ghosts = (src2->getGhostCellWidth()).max();

#if !defined(NDEBUG)
            if (V_ghosts != (src2->getGhostCellWidth()).min())
            {
                TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                           << "  src2 does not have uniform ghost cell widths" << std::endl);
            }

            if (patch_box != src2->getBox())
            {
                TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                           << "  dst, alpha, src1, and src2 must all live on the same patch" << std::endl);
            }
#endif
            MULTIPLY_ADD2_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, beta, V, V_ghosts, data_box.lower(0),
                             data_box.upper(0), data_box.lower(1), data_box.upper(1)
#if (NDIM == 3)
                                                                       ,
                             data_box.lower(2), data_box.upper(2)
#endif
                                 );
        }
    }
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseMultiply(Pointer<SideData<NDIM, double> > dst,
                                     const Pointer<SideData<NDIM, double> > alpha,
                                     const Pointer<SideData<NDIM, double> > src1,
                                     const Pointer<SideData<NDIM, double> > beta,
                                     const Pointer<SideData<NDIM, double> > src2,
                                     const Pointer<Patch<NDIM> > patch,
                                     const int i,
                                     const int j,
                                     const int k,
                                     const int l,
                                     const int m) const
{
    if (!src2)
    {
        const double zero = 0.0;
        pointwiseMultiply(dst, alpha, src1, zero, src2, patch, i, j, k, l);
        return;
    }

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        double* const D = dst->getPointer(axis, i);
        const int D_ghosts = (dst->getGhostCellWidth()).max();

        const double* const U = src1->getPointer(axis, j);
        const int U_ghosts = (src1->getGhostCellWidth()).max();

        const double* const V = src2->getPointer(axis, k);
        const int V_ghosts = (src2->getGhostCellWidth()).max();

        const double* const A = alpha->getPointer(axis, l);
        const int A_ghosts = (alpha->getGhostCellWidth()).max();

        const double* const B = beta->getPointer(axis, m);
        const int B_ghosts = (beta->getGhostCellWidth()).max();

        const Box<NDIM>& patch_box = patch->getBox();
        const Box<NDIM> data_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);

#if !defined(NDEBUG)
        if (D_ghosts != (dst->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst does not have uniform ghost cell widths" << std::endl);
        }

        if (U_ghosts != (src1->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src1 does not have uniform ghost cell widths" << std::endl);
        }

        if (V_ghosts != (src2->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  src2 does not have uniform ghost cell widths" << std::endl);
        }

        if (A_ghosts != (alpha->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  alpha does not have uniform ghost cell widths" << std::endl);
        }

        if (B_ghosts != (beta->getGhostCellWidth()).min())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  beta does not have uniform ghost cell widths" << std::endl);
        }

        if (patch_box != dst->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != src1->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != src2->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != alpha->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
        }

        if (patch_box != beta->getBox())
        {
            TBOX_ERROR("PatchMathOps::pointwiseMultiply():\n"
                       << "  dst, alpha, src1, beta, and src2 must all live on the same patch" << std::endl);
        }
#endif

        MULTIPLY_ADD3_FC(D, D_ghosts, A, A_ghosts, U, U_ghosts, B, B_ghosts, V, V_ghosts, data_box.lower(0),
                         data_box.upper(0), data_box.lower(1), data_box.upper(1)
#if (NDIM == 3)
                                                                   ,
                         data_box.lower(2), data_box.upper(2)
#endif
                             );
    }
    return;
} // pointwiseMultiply

void PatchMathOps::pointwiseL1Norm(Pointer<CellData<NDIM, double> > dst,
                                   const Pointer<CellData<NDIM, double> > src,
                                   const Pointer<Patch<NDIM> > patch) const
{
    double* const U = dst->getPointer();
    const int U_ghosts = (dst->getGhostCellWidth()).max();

    const double* const V = src->getPointer();
    const int V_ghosts = (src->getGhostCellWidth()).max();
    const int V_depth = src->getDepth();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (U_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL1Norm():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    const int U_depth = dst->getDepth();
    if (U_depth != 1)
    {
        TBOX_ERROR("PatchMathOps::pointwiseL1Norm():\n"
                   << "  dst does not have depth == 1" << std::endl);
    }

    if (V_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL1Norm():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL1Norm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL1Norm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    PW_L1_NORM_FC(U, U_ghosts, V, V_ghosts, V_depth, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                  patch_box.upper(1)
#if (NDIM == 3)
                      ,
                  patch_box.lower(2), patch_box.upper(2)
#endif
                      );
    return;
} // pointwiseL1Norm

void PatchMathOps::pointwiseL2Norm(Pointer<CellData<NDIM, double> > dst,
                                   const Pointer<CellData<NDIM, double> > src,
                                   const Pointer<Patch<NDIM> > patch) const
{
    double* const U = dst->getPointer();
    const int U_ghosts = (dst->getGhostCellWidth()).max();

    const double* const V = src->getPointer();
    const int V_ghosts = (src->getGhostCellWidth()).max();
    const int V_depth = src->getDepth();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (U_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL2Norm():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    const int U_depth = dst->getDepth();
    if (U_depth != 1)
    {
        TBOX_ERROR("PatchMathOps::pointwiseL2Norm():\n"
                   << "  dst does not have depth == 1" << std::endl);
    }

    if (V_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL2Norm():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL2Norm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL2Norm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    PW_L2_NORM_FC(U, U_ghosts, V, V_ghosts, V_depth, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                  patch_box.upper(1)
#if (NDIM == 3)
                      ,
                  patch_box.lower(2), patch_box.upper(2)
#endif
                      );
    return;
} // pointwiseL2Norm

void PatchMathOps::pointwiseMaxNorm(Pointer<CellData<NDIM, double> > dst,
                                    const Pointer<CellData<NDIM, double> > src,
                                    const Pointer<Patch<NDIM> > patch) const
{
    double* const U = dst->getPointer();
    const int U_ghosts = (dst->getGhostCellWidth()).max();

    const double* const V = src->getPointer();
    const int V_ghosts = (src->getGhostCellWidth()).max();
    const int V_depth = src->getDepth();

    const Box<NDIM>& patch_box = patch->getBox();

#if !defined(NDEBUG)
    if (U_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMaxNorm():\n"
                   << "  dst does not have uniform ghost cell widths" << std::endl);
    }

    const int U_depth = dst->getDepth();
    if (U_depth != 1)
    {
        TBOX_ERROR("PatchMathOps::pointwiseMaxNorm():\n"
                   << "  dst does not have depth == 1" << std::endl);
    }

    if (V_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMaxNorm():\n"
                   << "  src does not have uniform ghost cell widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMaxNorm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMaxNorm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    PW_MAX_NORM_FC(U, U_ghosts, V, V_ghosts, V_depth, patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                   patch_box.upper(1)
#if (NDIM == 3)
                       ,
                   patch_box.lower(2), patch_box.upper(2)
#endif
                       );
    return;
} // pointwiseMaxNorm

void PatchMathOps::pointwiseL1Norm(Pointer<NodeData<NDIM, double> > dst,
                                   const Pointer<NodeData<NDIM, double> > src,
                                   const Pointer<Patch<NDIM> > patch) const
{
    double* const U = dst->getPointer();
    const int U_ghosts = (dst->getGhostCellWidth()).max();

    const double* const V = src->getPointer();
    const int V_ghosts = (src->getGhostCellWidth()).max();
    const int V_depth = src->getDepth();

    const Box<NDIM>& patch_box = patch->getBox();
    const Box<NDIM> data_box = NodeGeometry<NDIM>::toNodeBox(patch_box);

#if !defined(NDEBUG)
    if (U_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL1Norm():\n"
                   << "  dst does not have uniform ghost node widths" << std::endl);
    }

    const int U_depth = dst->getDepth();
    if (U_depth != 1)
    {
        TBOX_ERROR("PatchMathOps::pointwiseL1Norm():\n"
                   << "  dst does not have depth == 1" << std::endl);
    }

    if (V_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL1Norm():\n"
                   << "  src does not have uniform ghost node widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL1Norm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL1Norm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    PW_L1_NORM_FC(U, U_ghosts, V, V_ghosts, V_depth, data_box.lower(0), data_box.upper(0), data_box.lower(1),
                  data_box.upper(1)
#if (NDIM == 3)
                      ,
                  data_box.lower(2), data_box.upper(2)
#endif
                      );
    return;
} // pointwiseL1Norm

void PatchMathOps::pointwiseL2Norm(Pointer<NodeData<NDIM, double> > dst,
                                   const Pointer<NodeData<NDIM, double> > src,
                                   const Pointer<Patch<NDIM> > patch) const
{
    double* const U = dst->getPointer();
    const int U_ghosts = (dst->getGhostCellWidth()).max();

    const double* const V = src->getPointer();
    const int V_ghosts = (src->getGhostCellWidth()).max();
    const int V_depth = src->getDepth();

    const Box<NDIM>& patch_box = patch->getBox();
    const Box<NDIM> data_box = NodeGeometry<NDIM>::toNodeBox(patch_box);

#if !defined(NDEBUG)
    if (U_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL2Norm():\n"
                   << "  dst does not have uniform ghost node widths" << std::endl);
    }

    const int U_depth = dst->getDepth();
    if (U_depth != 1)
    {
        TBOX_ERROR("PatchMathOps::pointwiseL2Norm():\n"
                   << "  dst does not have depth == 1" << std::endl);
    }

    if (V_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL2Norm():\n"
                   << "  src does not have uniform ghost node widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL2Norm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseL2Norm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    PW_L2_NORM_FC(U, U_ghosts, V, V_ghosts, V_depth, data_box.lower(0), data_box.upper(0), data_box.lower(1),
                  data_box.upper(1)
#if (NDIM == 3)
                      ,
                  data_box.lower(2), data_box.upper(2)
#endif
                      );
    return;
} // pointwiseL2Norm

void PatchMathOps::pointwiseMaxNorm(Pointer<NodeData<NDIM, double> > dst,
                                    const Pointer<NodeData<NDIM, double> > src,
                                    const Pointer<Patch<NDIM> > patch) const
{
    double* const U = dst->getPointer();
    const int U_ghosts = (dst->getGhostCellWidth()).max();

    const double* const V = src->getPointer();
    const int V_ghosts = (src->getGhostCellWidth()).max();
    const int V_depth = src->getDepth();

    const Box<NDIM>& patch_box = patch->getBox();
    const Box<NDIM> data_box = NodeGeometry<NDIM>::toNodeBox(patch_box);

#if !defined(NDEBUG)
    if (U_ghosts != (dst->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMaxNorm():\n"
                   << "  dst does not have uniform ghost node widths" << std::endl);
    }

    const int U_depth = dst->getDepth();
    if (U_depth != 1)
    {
        TBOX_ERROR("PatchMathOps::pointwiseMaxNorm():\n"
                   << "  dst does not have depth == 1" << std::endl);
    }

    if (V_ghosts != (src->getGhostCellWidth()).min())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMaxNorm():\n"
                   << "  src does not have uniform ghost node widths" << std::endl);
    }

    if (patch_box != dst->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMaxNorm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }

    if (patch_box != src->getBox())
    {
        TBOX_ERROR("PatchMathOps::pointwiseMaxNorm():\n"
                   << "  dst and src must live on the same patch" << std::endl);
    }
#endif

    PW_MAX_NORM_FC(U, U_ghosts, V, V_ghosts, V_depth, data_box.lower(0), data_box.upper(0), data_box.lower(1),
                   data_box.upper(1)
#if (NDIM == 3)
                       ,
                   data_box.lower(2), data_box.upper(2)
#endif
                       );
    return;
} // pointwiseMaxNorm

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
