# SGS Model API Sketch for `INSStaggeredHierarchyIntegrator`

This note sketches a more flexible design for subgrid-scale (SGS) turbulence
models in the staggered-grid incompressible Navier-Stokes solver. The main goal
is to separate the pieces that are shared by many models from the pieces that
actually distinguish one closure from another.

## Motivation

The current Smagorinsky implementation combines several different tasks in a
single class:

1. Fill ghost cells for the staggered velocity field
2. Reconstruct velocity-gradient information needed by the closure
3. Apply a constitutive closure to compute an SGS stress
4. Take the divergence of that stress and add it to the momentum equation

This works for one model, but it is not the most natural organization for a
larger family of closures. Most algebraic LES models differ primarily in how
they use velocity and velocity-gradient information, not in how ghost filling or
stress-divergence application is done.

## Terminology

It is reasonable to describe models as either `static` or `dynamic`:

- `static`: model coefficients are prescribed a priori
- `dynamic`: model coefficients are computed from the resolved flow, usually
  through test filtering and one or more averaging procedures

This terminology is common in the LES literature, but it should probably be
treated as a secondary property of a model rather than the primary software
split.

For software purposes, a more useful top-level classification is:

- algebraic eddy-viscosity models
- explicit stress models
- transport-equation SGS models

The first two fit naturally into the current explicit-force framework. The third
class will likely require additional state variables and should probably be
treated separately when we get there.

## Proposed Decomposition

The recommended split is:

1. Kinematics
2. Constitutive model
3. Stress application

### 1. Kinematics

This layer extracts resolved-flow quantities from the original staggered
velocity field. It is responsible for discrete finite-difference evaluation of
the quantities needed by the closure at their native locations, such as

- `grad(U)`
- strain-rate tensor `S`
- rotation tensor `Omega`
- invariants such as `|S|`
- optional filtered fields needed by dynamic or similarity models

`kinematics` is a good name for this layer: it describes quantities derived from
motion without yet imposing a turbulence closure.

### 2. Constitutive Model

This layer is the actual SGS closure. It consumes kinematic quantities and
produces either

- a dynamic eddy viscosity `mu_t`
- a staggered SGS stress tensor `tau_sgs`
- or both

Examples:

- Smagorinsky: compute `mu_t` from `|S|`
- WALE: compute `mu_t` from velocity-gradient combinations
- Vreman: compute `mu_t` from gradient invariants
- scale-similarity or mixed models: compute `tau_sgs` directly
- dynamic Smagorinsky: compute a coefficient from test-filtered quantities and
  then form `mu_t` or `tau_sgs`

For IBAMR, it is preferable to formulate the constitutive law in terms of the
dynamic viscosity `mu_t` rather than the kinematic viscosity `nu_t`. This fits
better with the rest of the library and naturally accommodates variable-density
extensions.

### 3. Stress Application

This layer converts a staggered SGS stress into the form needed by the momentum
equation, i.e.

`f_sgs = div(tau_sgs)`

on the staggered side-centered grid.

This part should be shared by most closures. It is largely a discrete operator
issue rather than a modeling issue.

## Recommended Common Product

For solver-facing SGS models, the common internal representation should probably
be a staggered stress tensor consistent with the MAC grid, not a uniformly
cell-centered symmetric tensor.

This is different from the statistics machinery. For turbulence statistics, it
is reasonable to pick a single analysis centering such as `CELL` or `NODE`
because the goal is postprocessing and visualization. For SGS closures, the
stress is part of the discrete operator, so its natural centering should follow
the staggered discretization.

That suggests the following layout:

- diagonal components such as `tau_xx`, `tau_yy`, `tau_zz` are cell-centered
- off-diagonal components such as `tau_xy`, `tau_xz`, `tau_yz` live on the
  staggered locations naturally associated with the discrete divergence

In practice this means the common product is a staggered stress container rather
than a single variable with one centering.

### Suggested Staggered Layout

In 2D:

- `tau_xx`, `tau_yy`: cell-centered
- `tau_xy`: node-centered

In 3D:

- `tau_xx`, `tau_yy`, `tau_zz`: cell-centered
- `tau_xy`, `tau_xz`, `tau_yz`: edge-centered

This is the tensor layout that is naturally dual to a side-centered momentum
equation:

- the `x`-momentum equation needs `d(tau_xx)/dx`, `d(tau_xy)/dy`, and
  `d(tau_xz)/dz`
- the `y`-momentum equation needs `d(tau_xy)/dx`, `d(tau_yy)/dy`, and
  `d(tau_yz)/dz`
- the `z`-momentum equation needs `d(tau_xz)/dx`, `d(tau_yz)/dy`, and
  `d(tau_zz)/dz`

Using the native staggered locations for these quantities should avoid
interpolating constitutive outputs, keep the discrete divergence more local, and
better match the boundary treatment and operator structure of the underlying
solver.

This also gives a common target for both eddy-viscosity and explicit-stress
models:

- eddy-viscosity models compute the staggered stress corresponding to
  `tau_sgs = 2 mu_t S`
- explicit-stress models write the staggered stress directly

This avoids over-specializing the interface around `mu_t`, which would be too
narrow for several important closures.

## Co-location Requirement

The SGS closure should be evaluated at the locations where its outputs are used.
This is especially important for algebraic models like WALE and Vreman, in
which `mu_t` depends nonlinearly on the resolved velocity-gradient tensor.

In particular:

- do not compute a single cell-centered `grad(U)` and reuse it everywhere
- do not interpolate resolved gradients from one centering to another
- do not compute `mu_t` at cell centers and interpolate it to node/edge
  locations
- do not compute a cell-centered SGS stress and then interpolate it back to the
  staggered stress layout

Instead:

- compute the resolved gradients directly at the locations where the
  corresponding stress component is defined
- use the finite-difference stencil appropriate to that target location on the
  staggered grid
- evaluate `mu_t` directly at those same locations
- form each `tau_ij = 2 mu_t S_ij` directly there

So, for example:

- the `mu_t` used in `tau_xx` should be computed from gradients co-located at
  cell centers
- the `mu_t` used in `tau_xy` should be computed from gradients co-located at
  the node/edge location where `tau_xy` is stored

This preserves the intended discrete model. In general, for nonlinear closures,
`compute then interpolate` is not equivalent to `interpolate the resolved data,
then compute`.

More strongly, `differentiate then interpolate` is not the same discrete
operation as `differentiate directly at the target location`. The kinematics
layer should therefore evaluate the needed derivatives from the original
staggered velocity components with the correct stencil for each target
location.

## Capability / Requirement Declaration

The SGS model should be able to declare what it needs. This lets the shared
infrastructure compute only the resolved quantities required by a specific
closure.

For example, a model may need:

- ghost-filled staggered `U`
- staggered diagonal strain components
- staggered off-diagonal strain components
- co-located invariants such as `|S|`
- a test filter
- global or plane averaging
- wall distance

A possible interface pattern is:

```cpp
enum class SGSKinematicsFlag
{
    U_GHOST_FILLED,
    STRAIN_RATE_DIAGONAL_CC,
    STRAIN_RATE_SHEAR_STAGGERED,
    STRAIN_RATE_MAGNITUDE_DIAGONAL_CC,
    STRAIN_RATE_MAGNITUDE_SHEAR_STAGGERED,
    TEST_FILTERED_U_CC,
    TEST_FILTERED_GRAD_U_CC,
    WALL_DISTANCE_CC
};

virtual std::set<SGSKinematicsFlag> getRequiredKinematics() const = 0;
```

Then a shared preparation layer can build only the requested fields before
calling the model. Those fields should be computed directly from the original
staggered velocity using the finite-difference stencil appropriate to the
requested location, not by interpolating gradients between centerings.

## Suggested Class Structure

One possible direction is:

### Top-level solver hook

Keep

- `INSTurbulenceModel`

as the object owned by `INSStaggeredHierarchyIntegrator`.

This class remains responsible for coordinating the full SGS contribution to the
momentum equation.

### Kinematics helper

Add a helper class, for example

- `INSSGSKinematics`

that

- fills ghost cells
- computes requested resolved quantities at the locations required by the
  staggered stress layout
- uses the finite-difference stencil appropriate to each target location
- manages scratch variables for those quantities

This helper should be reusable by multiple SGS closures.

### Constitutive base class

Add a base class for closures that produce a staggered SGS stress, for example

- `INSSGSStressModel`

with a method conceptually like

```cpp
virtual void computeSGSStress(
    const StaggeredSGSStressData& tau_data,
    const SGSKinematicsContext& context,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
    double data_time,
    const StokesSpecifications& problem_coefs) = 0;
```

where:

- `SGSKinematicsContext` provides access to the requested resolved quantities
- `StaggeredSGSStressData` bundles the cell-centered diagonal stresses and the
  node- or edge-centered shear stresses

Derived classes could include:

- `INSSmagorinskySGSStressModel`
- `INSWALESGSStressModel`
- `INSVremanSGSStressModel`
- `INSDynamicSmagorinskySGSStressModel`

### Shared stress-divergence wrapper

Then the top-level turbulence model can reuse one common implementation for:

- allocating `tau_sgs`
- asking the constitutive model for `tau_sgs`
- computing `div(tau_sgs)` on sides
- returning the explicit forcing to the integrator

## Static vs Dynamic in This Design

In this structure, `static` versus `dynamic` is not a separate top-level class
family. Instead, it is a property of a particular constitutive model.

For example:

- `Smagorinsky`: static algebraic eddy-viscosity model
- `DynamicSmagorinsky`: dynamic algebraic eddy-viscosity model

Both are stress models. The dynamic version simply declares additional
requirements such as a test filter and averaging machinery.

## Near-Term Refactor Path

A low-risk path from the current code would be:

1. Extract common ghost-fill and stress-divergence logic from the current
   Smagorinsky model
2. Introduce a small kinematics helper for co-located staggered
   strain/gradient quantities computed directly with native-location
   finite-difference stencils
3. Introduce a staggered-stress container with cell-centered diagonal and
   node-/edge-centered shear components
4. Refactor the current Smagorinsky class into a constitutive stress model that
   evaluates `mu_t` and `tau_sgs` directly at those native locations
5. Add one more algebraic model, such as WALE, to validate the design
6. Add a dynamic Smagorinsky model only after the test-filter and averaging
   interfaces are in place

This should give us a path that is flexible enough for the next several LES
closures without prematurely designing around the hardest transport-equation
models.
