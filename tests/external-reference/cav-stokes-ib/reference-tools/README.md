# External CAV Stokes-IB Reference Generation (Octave Baseline)

This directory stores scripts for generating local reference artifacts from the
MATLAB/Octave baseline model.

## Recommended workflow

1. Build IBAMR tests.
2. Run:

```bash
tests/external-reference/cav-stokes-ib/reference-tools/generate_single_level_case_01.sh \
  <build_dir> <path_to_implicit_ib_coupling_aware_vanka>
```

For all external-reference inputs in this directory (single-level, multilevel_l2,
multilevel_l3), run:

```bash
tests/external-reference/cav-stokes-ib/reference-tools/generate_all_cases.sh \
  <build_dir> <path_to_implicit_ib_coupling_aware_vanka>
```

The script performs two steps:
1. Runs the IBAMR parity executable once in `write_reference = TRUE` mode to export `eulerian_dof_map.ref`.
2. Runs Octave to generate independent MATLAB-baseline artifacts (`X/A/F/S/J/SAJ/A00`) reordered into IBAMR Eulerian DOF numbering using that map.

The helper script rewrites a temporary input file so `reference_dir` is an
absolute build-tree path.

By default the script writes local artifacts under:

- `<build_dir>/tests/external-reference/cav-stokes-ib/reference-generated/single_level_case_01/`

Do not commit generated `.ref` artifacts.
