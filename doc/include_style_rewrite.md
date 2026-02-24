# Include Style Rewrite Script

Use this script to normalize include directives between quoted and angle-bracket style:

- `#include "foo/bar.h"` -> `#include <foo/bar.h>`
- `#include <MyLocalHeader.h>` -> `#include "MyLocalHeader.h"` when `MyLocalHeader.h` exists in the same directory as the source file.

The default exception keeps `tests.h` includes quoted (e.g. `"../tests.h"`), and same-directory private headers/sources are also kept quoted automatically.
The script also skips `ibtk/include/samrai_compatibility/legacy/*` wrappers.

## Script

- `/Users/boyceg/code/IBAMR/scripts/maintenance/rewrite_include_quotes_to_angle_brackets`

## Dry-run

```bash
/Users/boyceg/code/IBAMR/scripts/maintenance/rewrite_include_quotes_to_angle_brackets \
  /Users/boyceg/code/IBAMR/examples \
  /Users/boyceg/code/IBAMR/include \
  /Users/boyceg/code/IBAMR/src \
  /Users/boyceg/code/IBAMR/ibtk/examples \
  /Users/boyceg/code/IBAMR/ibtk/include \
  /Users/boyceg/code/IBAMR/ibtk/src \
  /Users/boyceg/code/IBAMR/tests
```

## Apply changes

```bash
/Users/boyceg/code/IBAMR/scripts/maintenance/rewrite_include_quotes_to_angle_brackets --write \
  /Users/boyceg/code/IBAMR/examples \
  /Users/boyceg/code/IBAMR/include \
  /Users/boyceg/code/IBAMR/src \
  /Users/boyceg/code/IBAMR/ibtk/examples \
  /Users/boyceg/code/IBAMR/ibtk/include \
  /Users/boyceg/code/IBAMR/ibtk/src \
  /Users/boyceg/code/IBAMR/tests
```

## Optional custom exception regex

If you need additional quoted exceptions:

```bash
/Users/boyceg/code/IBAMR/scripts/maintenance/rewrite_include_quotes_to_angle_brackets --write \
  --keep-quoted-regex '(^|/)(tests\\.h|my_local_header\\.h)$' \
  /Users/boyceg/code/IBAMR/tests
```
