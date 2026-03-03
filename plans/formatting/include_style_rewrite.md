# Include Style Rewrite Script

This document describes `scripts/maintenance/rewrite_include_quotes_to_angle_brackets.py`.
Use that script to normalize include directives between quoted and angle-bracket style:

- `#include "foo/bar.h"` -> `#include <foo/bar.h>`
- `#include <MyLocalHeader.h>` -> `#include "MyLocalHeader.h"` when `MyLocalHeader.h` exists in the same directory as the source file.

The default exception keeps `tests.h` includes quoted (e.g. `"../tests.h"`), and same-directory private headers/sources are also kept quoted automatically.

## Script

- `scripts/maintenance/rewrite_include_quotes_to_angle_brackets.py`

## Canonical entrypoints

From an out-of-tree CMake build directory:

```bash
cd /path/to/build-dir
make indent
```

or from anywhere:

```bash
cmake --build /path/to/build-dir --target indent
```

This applies include-style rewrite to recently changed files before clang-format and other normalization steps.

From an out-of-tree CMake build directory:

```bash
cd /path/to/build-dir
make indent-all
```

or from anywhere:

```bash
cmake --build /path/to/build-dir --target indent-all
```

This applies include-style rewrite to the full tree before clang-format and other normalization steps.

## Dry-run

Run from the repository root:

```bash
./scripts/maintenance/rewrite_include_quotes_to_angle_brackets.py \
  examples \
  include \
  src \
  ibtk/examples \
  ibtk/include \
  ibtk/src \
  tests
```

## Apply changes

Run from the repository root:

```bash
./scripts/maintenance/rewrite_include_quotes_to_angle_brackets.py --write \
  examples \
  include \
  src \
  ibtk/examples \
  ibtk/include \
  ibtk/src \
  tests
```

## Optional custom exception regex

If you need additional quoted exceptions:
Run from the repository root:

```bash
./scripts/maintenance/rewrite_include_quotes_to_angle_brackets.py --write \
  --keep-quoted-regex '(^|/)(tests\\.h|my_local_header\\.h)$' \
  tests
```
