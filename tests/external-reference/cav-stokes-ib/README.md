## External Reference Verification: CAV Stokes-IB

This directory contains developer-only tests that verify IBAMR behavior against
an external MATLAB/Octave reference implementation.

These tests are not part of routine CI.

### CI vs. dev-only

- CI keeps only self-contained tests in `tests/IB` and `tests/navier_stokes`.
- Files in this directory require additional infrastructure and are run
  manually.

### Run the dev-only external-reference test

```bash
cmake --build <build_dir> -j10 --target run-external-reference-cav-stokes-ib
```

### Run the dev-only external-reference test via attest

```bash
cmake --build <build_dir> -j10 --target run-external-reference-cav-stokes-ib-attest
```

This runs all attest cases with prefix:

- `external-reference/cav-stokes-ib/implicit_stokes_ib_reference_single_level_01`

If local references are missing, generate them first:

```bash
tests/external-reference/cav-stokes-ib/reference-tools/generate_single_level_case_01.sh \
  <build_dir> <path_to_implicit_ib_coupling_aware_vanka>
```

To regenerate all external-reference cases (single-level + multilevel):

```bash
IBAMR_EXTERNAL_VANKA_REPO=<path_to_implicit_ib_coupling_aware_vanka> \
cmake --build <build_dir> -j10 --target generate-external-reference-cav-stokes-ib
```

### Optional reference generation workflow

Use the scripts in `reference-tools/` to generate local references with Octave
and the external MATLAB baseline repository.

The run/generation scripts rewrite temporary input files so `reference_dir`
points to the build-tree external-reference directory.

Generated `.ref` files should remain local and must not be committed.
By default these external-reference inputs write/read references under the
build directory (`build/tests/external-reference/cav-stokes-ib/reference-generated/...`).
