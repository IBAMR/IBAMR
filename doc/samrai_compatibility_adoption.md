# Adopting the SAMRAI Compatibility Layer in Another Codebase

This note describes two migration paths: scripted conversion and manual conversion.

## Scripted Conversion

Use the rewrite script for broad, mechanical updates across many files.

1. Build a deterministic file list (for example, headers and C++ sources).
2. Run a dry-run first.
3. Run with `--write` plus safety checks.
4. Build and fix edge cases manually.

Small example:

```bash
# 1) Candidate file list
find src include -type f \( -name '*.h' -o -name '*.cpp' \) > /tmp/samrai_targets.txt

# 2) Dry-run
python3 scripts/maintenance/rewrite_samrai_compatibility_file.py $(cat /tmp/samrai_targets.txt)

# 3) Apply with safety checks
python3 scripts/maintenance/rewrite_samrai_compatibility_file.py \
  --write \
  --fail-on-string-change \
  --fail-on-samrai-forward-decls \
  $(cat /tmp/samrai_targets.txt)
```

Typical script output change in one file:

```cpp
// Before
#include <CellData.h>
SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > u;

// After
#include "SAMRAICellData.h"
#include "SAMRAIPointer.h"
SAMRAIPointer<SAMRAICellData<double> > u;
```

Safety guarantees enforced by the script:

1. String literals are not allowed to change (`--fail-on-string-change`).
2. SAMRAI forward declarations are rewritten away and must not remain (`--fail-on-samrai-forward-decls`).
3. Pointer kind is preserved:
`Pointer<T> -> SAMRAIPointer<T>`, `ConstPointer<T> -> SAMRAIConstPointer<T>`.

## Manual Adoption Process

Use this for small migrations or follow-up cleanup after scripted changes.

1. Ensure include paths contain `ibtk/include/ibtk` and `ibtk/include/samrai_compatibility`.
2. Add `"ibtk/samrai_compatibility_names.h"`.
3. Replace direct SAMRAI includes with specific `SAMRAI*.h` compatibility includes.
4. Replace `SAMRAI::...::Type` (or unscoped legacy names) with `SAMRAIType` aliases.
5. Build immediately and fix collisions one file at a time.

Small example:

```cpp
// Before
#include <IntVector.h>
#include <Box.h>
#include <tbox/Pointer.h>
#include <CellData.h>

SAMRAI::hier::IntVector<NDIM> ghosts(1);
SAMRAI::hier::Box<NDIM> box = patch->getBox();
SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > data =
    patch->getPatchData(data_idx);

// After
#include "ibtk/samrai_compatibility_names.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIBox.h"
#include "SAMRAIPointer.h"
#include "SAMRAICellData.h"

SAMRAIIntVector ghosts(1);
SAMRAIBox box = patch->getBox();
SAMRAIPointer<SAMRAICellData<double> > data = patch->getPatchData(data_idx);
```
