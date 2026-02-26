# Adopting the SAMRAI Compatibility Layer in Another Codebase

This note describes two migration paths: scripted conversion and manual conversion.

## Scripted Conversion

Use the rewrite scripts for broad, mechanical updates across many files.

From an out-of-tree CMake build directory:

```bash
cd /path/to/build-dir
make rewrite-samrai-compatibility
```

or from anywhere:

```bash
cmake --build /path/to/build-dir --target rewrite-samrai-compatibility
```

Apply only to recently changed files.

From an out-of-tree CMake build directory:

```bash
cd /path/to/build-dir
make rewrite-samrai-compatibility-all
```

or from anywhere:

```bash
cmake --build /path/to/build-dir --target rewrite-samrai-compatibility-all
```

Apply to the full tree (`include`, `src`, `ibtk/include`, `ibtk/src`, `examples`, `ibtk/examples`, `tests`).

1. Iteratively apply SAMRAI alias/type rewrites until convergence.
2. Iteratively apply include quote-to-angle rewrites until convergence.
3. Require pre-format dry-run convergence (`updated=0` and `changed=0`).
4. Run formatting (`make indent-all` in the build directory).
5. Re-run both scripts in dry-run mode for post-format reporting (informational).
6. Build and fix remaining one-off compile issues.

Small example:

```bash
# 1) Candidate file list
find src include -type f \( -name '*.h' -o -name '*.cpp' \) > /tmp/samrai_targets.txt

# 2) Alias rewrite dry-run
python3 scripts/maintenance/rewrite_samrai_compatibility_file.py $(cat /tmp/samrai_targets.txt)

# 3) Alias rewrite apply with strict checks
python3 scripts/maintenance/rewrite_samrai_compatibility_file.py \
  --write \
  --fail-on-string-change \
  --fail-on-samrai-forward-decls \
  $(cat /tmp/samrai_targets.txt)

# 4) Include quote-to-angle dry-run
python3 scripts/maintenance/rewrite_include_quotes_to_angle_brackets.py $(cat /tmp/samrai_targets.txt)

# 5) Include quote-to-angle apply
python3 scripts/maintenance/rewrite_include_quotes_to_angle_brackets.py \
  --write \
  $(cat /tmp/samrai_targets.txt)
```

Notes:

1. The adoption runner script performs both rewrites iteratively to convergence before formatting.
2. The required convergence gate is pre-format dry-run (`updated=0` and `changed=0`).
3. Post-format dry-runs are reported but non-fatal, since `indent-all` may alter include grouping/spacing.
4. The scripts do not build; run the appropriate CMake/autotools build gates afterward.
5. Commit-ready source must match repository formatting rules; always run `make indent-all` before committing.

Typical conversion in one file:

```cpp
// Before
#include <CellData.h>
SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> u;

// After
#include <SAMRAICellData.h>
#include <SAMRAIPointer.h>
SAMRAIPointer<SAMRAICellData<double>> u;
```

Strict safety checks are available when you pass specific flags:

1. `--fail-on-string-change`: fail if any string literal changes.
2. `--fail-on-samrai-forward-decls`: fail if SAMRAI forward declarations remain.
3. Pointer kind preservation: `Pointer<T> -> SAMRAIPointer<T>`, `ConstPointer<T> -> SAMRAIConstPointer<T>`.

## Manual Adoption Process

Use this for small migrations or follow-up cleanup after scripted changes.

1. Ensure include paths contain `ibtk/include/ibtk` and `ibtk/include/samrai_compatibility`.
2. Add `<ibtk/samrai_compatibility_names.h>`.
3. Replace direct SAMRAI includes with specific `SAMRAI*.h` compatibility includes.
4. Replace `SAMRAI::...::Type` (or unscoped legacy names) with `SAMRAIType` aliases.
5. Remove SAMRAI forward declarations and include the needed compatibility headers instead.
6. Build immediately and fix collisions one file at a time.

Small example:

```cpp
// Before
#include <IntVector.h>
#include <Box.h>
#include <tbox/Pointer.h>
#include <CellData.h>

SAMRAI::hier::IntVector<NDIM> ghosts(1);
SAMRAI::hier::Box<NDIM> box = patch->getBox();
SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> data =
    patch->getPatchData(data_idx);

// After
#include <ibtk/samrai_compatibility_names.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIBox.h>
#include <SAMRAIPointer.h>
#include <SAMRAICellData.h>

SAMRAIIntVector ghosts(1);
SAMRAIBox box = patch->getBox();
SAMRAIPointer<SAMRAICellData<double>> data = patch->getPatchData(data_idx);
```
