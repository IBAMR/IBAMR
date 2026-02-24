#!/usr/bin/env python3

# ---------------------------------------------------------------------
#
# Copyright (c) 2026 by the IBAMR developers
# All rights reserved.
#
# This file is part of IBAMR.
#
# IBAMR is free software and is distributed under the 3-clause BSD
# license. The full text of the license can be found in the file
# COPYRIGHT at the top level directory of IBAMR.
#
# ---------------------------------------------------------------------

"""
Generate legacy SAMRAI include wrappers under
ibtk/include/samrai_compatibility/legacy.

Each wrapper preserves an old include name (e.g. CellVariable.h or
tbox/Pointer.h) and forwards to the corresponding compatibility header
SAMRAI*.h.
"""

from __future__ import annotations

import re
from pathlib import Path


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    alias_root = repo / "ibtk/include/samrai_compatibility"
    legacy_root = alias_root / "legacy"
    legacy_root.mkdir(parents=True, exist_ok=True)

    module_map: dict[str, str] = {}
    flat_candidates: dict[str, set[str]] = {}

    for header in alias_root.rglob("SAMRAI*.h"):
        rel = header.relative_to(alias_root)
        if len(rel.parts) < 2:
            continue
        module = rel.parts[0]
        stem = header.stem
        if not stem.startswith("SAMRAI"):
            continue
        old_name = f"{stem[len('SAMRAI'):]}.h"
        module_map[f"{module}/{old_name}"] = f"{stem}.h"
        flat_candidates.setdefault(old_name, set()).add(f"{stem}.h")

    flat_map = {k: next(iter(v)) for k, v in flat_candidates.items() if len(v) == 1}

    # Generate wrappers for all known aliases so every SAMRAI*.h shim has
    # a corresponding legacy include path.
    needed: set[str] = set(module_map.keys()) | set(flat_map.keys())

    count = 0
    for include_name in sorted(needed):
        shim_name = module_map.get(include_name) or flat_map.get(include_name)
        if not shim_name:
            continue

        out = legacy_root / include_name
        out.parent.mkdir(parents=True, exist_ok=True)

        guard = "included_IBTK_samrai_compatibility_legacy_" + re.sub(r"[^A-Za-z0-9]+", "_", include_name)
        text = f"""// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef {guard}
#define {guard}

#include <samrai_compatibility/samrai_compatibility_environment.h>
#include <{shim_name}>

#endif // #ifndef {guard}
"""
        out.write_text(text)
        count += 1

    print(f"generated {count} wrapper header(s) in {legacy_root}")


if __name__ == "__main__":
    main()
