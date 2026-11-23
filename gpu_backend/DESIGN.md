# GPU Fluid Prototype Design Note (For Reference — Not For Merge)

This branch contains our working 3D CUDA prototype for the conservative variable-coefficient staggered Navier-Stokes fluid solver, shared for architectural reference in support of IBAMR's SAMRAI 4 / IBAMR 1.0 GPU design.

## 1. Layout and Touchpoints

- **Modified IBAMR files (minimal hooks):**
  - `include/ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h`: adds the 3D bridge pointer and cached per-regrid flags (`d_gpu_bridge`, `d_gpu_hierarchy_dirty`, `d_gpu_brinkman_dirty`).
  - `src/navier_stokes/INSVCStaggeredConservativeHierarchyIntegrator.cpp`: adds opt-in eligibility check `canUseGpuCandidate2`, preprocess offload, Stokes solve offload, and safe CPU fallback path.
  - `src/CMakeLists.txt`: links `${CMAKE_DL_LIBS}` into `IBAMR2d` and `IBAMR3d`.
- **New bridge files in `src/navier_stokes/`:**
  - `INSVCStaggeredConservativeGPUBridge.h`: header-only `dlopen`/`dlsym` loader and MPI rank-0 gather/scatter bridge.
  - `ibamr_c2_gpu.h`: shared pure C ABI declaration.
- **Standalone CUDA library in `gpu_backend/`:**
  - `ibamr_c2_gpu.h`: pure C ABI (`C2GpuStatus`, `C2GridConfig`, `C2PatchBox`, `C2Telemetry`, and 9 `c2_gpu_*` functions).
  - `c2_gpu_kernels.cuh`: `__device__`/`__global__` inline kernels (no `throw`, no `virtual` calls).
  - `c2_gpu_engine.cu`: CUDA context management, 2-level ghost sync, SSPRK3 advection, and Krylov + geometric multigrid Stokes solver.
  - `test_c2_gpu_lib.cpp` & `build.sh`: self-contained library build and synthetic unit/capacity test.

## 2. Pure C Interface (`ibamr_c2_gpu.h`)

The CUDA engine is completely decoupled from SAMRAI/PETSc headers behind plain POD structs (`C2GridConfig`, `C2PatchBox`, `C2Telemetry`), `C2GpuStatus` integer return codes (`C2_OK`, `C2_ERR_INVALID`, `C2_ERR_CAPACITY`, `C2_ERR_CUDA`, `C2_ERR_NONFINITE`), and 9 C functions:

- Lifecycle & AMR layout: `c2_gpu_create`, `c2_gpu_destroy`, `c2_gpu_regrid`
- Static/state uploads: `c2_gpu_upload_brinkman`, `c2_gpu_upload_rho`, `c2_gpu_upload_state`
- Time-step compute: `c2_gpu_preprocess` (SSPRK3 advection + RHS assembly), `c2_gpu_integrate_solve` (outer Krylov + MG Stokes solve)
- Diagnostics: `c2_gpu_get_telemetry`

## 3. Bridge into `INSVCStaggeredConservativeHierarchyIntegrator`

When `use_gpu_fluid_backend = TRUE` is set in the integrator input database (default `FALSE`), `INSVCStaggeredConservativeGPUBridge` loads `libibamr_c2_gpu.so` on MPI rank 0 via `dlopen` (from `gpu_fluid_backend_db { library_path = ... }` or the executable directory via `dladdr`). At each time step, `canUseGpuCandidate2` checks structural preconditions (3D, at most 2 AMR levels with refinement ratio 2, `CUI` limiters, `SSPRK3` density stepping). If initialization, capacity, or solve returns a non-`C2_OK` status and `fallback_to_cpu = TRUE`, the integrator cleanly restores CPU preprocessing/Stokes state; with `fallback_to_cpu = FALSE`, it fails closed.

## 4. Host / Device Data Movement Per Step

- **On AMR regrid (`c2_gpu_regrid`, rank 0 gather once per regrid):** uploads the `L1` patch box list (`hi`/`lo` indices, slot offsets) and precomputes the coarse-to-fine owner lookup table on device.
- **Step entry (`c2_gpu_upload_brinkman` when dirty, `c2_gpu_upload_state`):** rank 0 gathers side-centered `(u, rho)` and cell-centered `p` across `L0` and `L1` and transfers packed arrays H2D; if the Eulerian state is uniform constant density without Brinkman changes, intermediate transfers are elided or kept compact.
- **SSPRK3 & RHS (`c2_gpu_preprocess`):** runs 3-stage SSPRK3 conservative momentum/density advection, coarse-fine restriction/prolongation, quadratic CF interpolation, and physical BCs entirely on GPU.
- **IB coupling & Stokes solve (`c2_gpu_integrate_solve`):** Lagrangian IBFE force evaluation, spreading, and velocity interpolation stay on CPU via existing SAMRAI/libMesh routines. Only the spread Eulerian force `f` is uploaded H2D; the matrix-free GPU Stokes solve computes updated `(u_new, p_new)` (and optional `rho_new`) on device and copies them D2H back to SAMRAI patch data.

## 5. Mixed-Precision FP32 Multigrid Inside FP64 Krylov

The Stokes solver keeps outer Krylov vectors, matrix-free residual evaluation (`A x`), and dot products in **FP64**, while optionally executing the geometric multigrid (GMG) V-cycle preconditioner in **FP32** (`enable_fp32_mg = TRUE`) to double V-cycle bandwidth and cut VRAM footprint. Each outer solve checks the true FP64 relative residual; if FP32 preconditioning stagnates or fails finite checks, it falls back to FP64 GMG (with optional `sticky_fp64_on_fb = TRUE` to hold FP64 MG for remaining steps). Both pure FP64 (`enable_fp32_mg = FALSE`) and FP32-MG modes are supported.

## 6. Current Limits

- 3D only, uniform `L0` plus at most one `L1` AMR level (`ratio = 2`).
- Single-GPU execution driven through MPI rank-0 gather/scatter (`MPI_Gatherv`/`MPI_Scatterv`).
- Channel-oriented physical boundary kernels (z-midplane symmetry + inlet/outlet pressure-velocity BCs).
- CUDA-only (`libibamr_c2_gpu.so` via `nvcc`) with static internal L1 slot capacity limits (returning `C2_ERR_CAPACITY` if exceeded).
- IBFE solid mechanics, IB spread/interp, and SAMRAI hierarchy regridding stay on CPU.

## 7. Measured Benchmark Results (8 CPU MPI Ranks + 8 GB NVIDIA GPU)

1. **3D channel-constriction FSI (224x50x25 base grid, 2 AMR levels, 5,120 steps, GPU FP32-MG vs matched CPU):**
   - **End-to-end speedup:** **`2.33x`** (`777.3 s` wall time vs `1807.8 s` CPU).
   - **Component breakdown:** fluid Stokes solve `694.4 s -> 138.5 s`; SSPRK3 advection `397.5 s -> 95.2 s`; peak GPU VRAM `1,342 MiB`; `0` CPU or FP64 fallbacks.
   - **Physics parity (129 samples):** max relative difference `< 0.4%` (`0.389%` x-drag force, `0.023%` max velocity, `0.050%` inlet pressure, `0.0073%` centroid).
2. **Changed-mesh reference-reset FSI leg (300 steps, GPU pure FP64 vs matched CPU):**
   - **Speedup:** **`2.51x`** (`53.5 s` vs `134.2 s` over 299 matched steps).
   - **FP64 parity:** max relative difference `3.6e-7` on cell centroid, `1.3e-6` on cell volume, `1.2e-4` on min quadrature Jacobian, `7.4e-4` on x-drag force, zero inversions, zero fallbacks.

## 8. Potential Reuse for SAMRAI 4 / IBAMR 1.0

- The device kernels (`c2_gpu_kernels.cuh`) already avoid C++ exceptions, virtual dispatch, and dynamic device allocation, using plain POD structs and inline device functions.
- Under SAMRAI 4 device patch allocation (`umpire` / Kokkos / CUDA memory spaces), the rank-0 MPI gather/scatter bridge in `INSVCStaggeredConservativeGPUBridge.h` can be eliminated in favor of direct per-patch device pointers while reusing the staggered SSPRK3, conservative CUI interpolation, matrix-free Stokes GMG V-cycle, and FP32-MG + FP64-Krylov pattern.
