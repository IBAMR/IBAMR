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
Rewrite a source/header file to use SAMRAI compatibility aliases.

This tool:
1) rewrites SAMRAI include directives to compatibility include directives
   (e.g. <CellVariable.h> -> <SAMRAICellVariable.h>),
2) rewrites SAMRAI type spellings to compatibility aliases
   (e.g. SAMRAI::hier::IntVector, hier::IntVector, IntVector -> SAMRAIIntVector),
3) optionally writes results in place.

The script is heuristic by design and may require manual follow-up for name
collisions (e.g. generic names like Index).
"""

from __future__ import annotations

import argparse
import dataclasses
import difflib
import pathlib
import re
import sys
from typing import Dict, List, Sequence, Set, Tuple


SAMRAI_MODULES = ("algs", "appu", "geom", "hier", "math", "mesh", "pdat", "solv", "tbox", "xfer")
SAMRAI_MODULES_PATTERN = "|".join(SAMRAI_MODULES)
USING_DECL_RE = re.compile(r"^\s*using\s+(SAMRAI[A-Za-z0-9_]+)\s*=\s*([^;]+);", re.MULTILINE)
INCLUDE_RE = re.compile(r'^(\s*#\s*include\s*)([<"])([^>"]+)([>"])(.*)$')
FORWARD_DECL_RE = re.compile(
    r"^\s*(?:template\s*<[^;{}]+>\s*)?(?:class|struct)\s+[A-Za-z_][A-Za-z0-9_]*\s*;\s*$"
)
USING_ALIAS_LINE_RE = re.compile(r"^\s*using\s+[A-Za-z_][A-Za-z0-9_]*\s*=\s*[^;]+;\s*$")
TEMPLATE_INSTANTIATION_RE = re.compile(r"^\s*template\s+(?:class|struct)\s+[^;]+;\s*$")
NON_CODE_RE = re.compile(
    r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',
    re.MULTILINE | re.DOTALL,
)
PROTECTED_CODE_TOKENS = (
    "IBTK::FACPreconditioner",
    "IBTK::FACPreconditionerStrategy",
    "IBMethod",
    "IBStrategy",
    "IBImplicitStrategy",
)
TYPE_TOKEN_RE = re.compile(
    r"(?:SAMRAI::)?(?:algs|appu|geom|hier|math|mesh|pdat|solv|tbox|xfer)::([A-Za-z_][A-Za-z0-9_]*)"
)
AMBIGUOUS_UNQUALIFIED_BASES = {
    "FACPreconditioner",
    "FACPreconditionerStrategy",
}


@dataclasses.dataclass
class AliasEntry:
    alias_name: str
    include_name: str
    include_keys: Set[str]
    base_types: Set[str]
    is_template_alias: bool


def _repo_root(script_path: pathlib.Path) -> pathlib.Path:
    return script_path.resolve().parents[2]


def _compat_dir(repo_root: pathlib.Path) -> pathlib.Path:
    return repo_root / "ibtk" / "include" / "samrai_compatibility"


def _load_alias_catalog(repo_root: pathlib.Path) -> List[AliasEntry]:
    compat_dir = _compat_dir(repo_root)
    entries: Dict[str, AliasEntry] = {}

    for header_path in compat_dir.rglob("SAMRAI*.h"):
        if "legacy" in header_path.parts:
            continue
        alias_name = header_path.stem
        rel_parent = header_path.parent.relative_to(compat_dir)
        include_name = f"{alias_name}.h"
        include_keys: Set[str] = set()
        base_types: Set[str] = set()
        text = header_path.read_text(encoding="utf-8")
        is_template_alias = False
        lines = text.splitlines()
        using_pat = re.compile(rf"^\s*using\s+{re.escape(alias_name)}\s*=")
        for i, line in enumerate(lines):
            if not using_pat.match(line):
                continue
            lo = max(0, i - 8)
            window = "\n".join(lines[lo:i])
            if re.search(r"template\s*<[^>]+>", window):
                is_template_alias = True
                break
        for line in text.splitlines():
            m = INCLUDE_RE.match(line)
            if not m:
                continue
            inc = m.group(3).strip()
            # Only map header-style includes. This prevents accidental
            # rewrites of standard library includes like <memory>.
            if not inc.endswith(".h"):
                continue
            if inc.startswith("ibtk/") or inc.startswith("samrai_compatibility/"):
                continue
            if inc == "SAMRAI_config.h":
                continue
            include_keys.add(inc)
            include_keys.add(pathlib.PurePosixPath(inc).name)

        for m in USING_DECL_RE.finditer(text):
            rhs = m.group(2)
            if "SAMRAI::" not in rhs:
                continue
            for t in TYPE_TOKEN_RE.findall(rhs):
                # Keep SAMRAI*-prefixed base names too (e.g. SAMRAIVectorReal):
                # these are valid SAMRAI types that still need alias rewrites,
                # typically to aliases like SAMRAISAMRAIVectorReal.
                base_types.add(t)

        if not base_types:
            continue

        # Synthetic include keys for SAMRAI4 style include paths.
        for base in base_types:
            include_keys.add(f"{base}.h")
            for module in SAMRAI_MODULES:
                include_keys.add(f"{module}/{base}.h")
                include_keys.add(f"SAMRAI/{module}/{base}.h")

        if alias_name not in entries:
            entries[alias_name] = AliasEntry(alias_name, include_name, include_keys, base_types, is_template_alias)
        else:
            entries[alias_name].include_keys.update(include_keys)
            entries[alias_name].base_types.update(base_types)
            entries[alias_name].is_template_alias = entries[alias_name].is_template_alias or is_template_alias

    # Longest alias names first for deterministic processing.
    return sorted(entries.values(), key=lambda e: len(e.alias_name), reverse=True)


def _replace_includes(
    lines: List[str], entries: Sequence[AliasEntry]
) -> Tuple[List[str], Dict[str, Set[str]], bool]:
    include_key_to_alias: Dict[str, Tuple[str, int]] = {}

    def _score_include_key(key: str, alias_name: str) -> int:
        key_stem = pathlib.PurePosixPath(key).name
        if key_stem.endswith(".h"):
            key_stem = key_stem[:-2]
        alias_stem = alias_name
        if alias_stem.startswith("SAMRAI"):
            alias_stem = alias_stem[len("SAMRAI") :]
        if key_stem == alias_stem:
            return 3
        if key_stem == alias_name:
            return 2
        return 1

    for entry in entries:
        for key in entry.include_keys:
            cand = (entry.include_name, _score_include_key(key, entry.alias_name))
            prev = include_key_to_alias.get(key)
            if prev is None or cand[1] > prev[1]:
                include_key_to_alias[key] = cand

    updated: List[str] = []
    changed = False
    file_include_keys: Set[str] = set()
    rewritten_aliases: Set[str] = set()

    for line in lines:
        m = INCLUDE_RE.match(line)
        if not m:
            updated.append(line)
            continue
        prefix, open_delim, inc, close_delim, suffix = m.groups()
        key = inc.strip()
        file_include_keys.add(key)
        file_include_keys.add(pathlib.PurePosixPath(key).name)
        if key.startswith("ibtk/"):
            updated.append(line)
            continue
        alias_inc_info = include_key_to_alias.get(key)
        alias_inc = alias_inc_info[0] if alias_inc_info is not None else None
        if alias_inc is None:
            alias_inc_info = include_key_to_alias.get(pathlib.PurePosixPath(key).name)
            alias_inc = alias_inc_info[0] if alias_inc_info is not None else None

        if alias_inc is not None and key != alias_inc:
            # Use quotes for compatibility headers by default.
            updated.append(f'{prefix}"{alias_inc}"{suffix}')
            rewritten_aliases.add(alias_inc[:-2])  # drop .h
            changed = True
        else:
            updated.append(line)

    return updated, {"file_include_keys": file_include_keys, "rewritten_aliases": rewritten_aliases}, changed


def _contains_qualified_type(text: str, base: str) -> bool:
    return bool(
        re.search(rf"\bSAMRAI::(?:[A-Za-z_][A-Za-z0-9_]*::)+{re.escape(base)}\b", text)
        or re.search(rf"\b(?:{SAMRAI_MODULES_PATTERN})::{re.escape(base)}\b", text)
    )


def _contains_ndim_template_use(text: str, base: str) -> bool:
    # Strong signal for SAMRAI data types that historically carry NDIM.
    # This avoids broad unqualified replacement while still catching common
    # forms like:
    #   Pointer<PatchHierarchy<NDIM> >
    #   SAMRAIPointer<CartesianGridGeometry<NDIM> >
    return bool(re.search(rf"\b{re.escape(base)}\s*<\s*NDIM\b", text))


def _replace_types(
    text: str, entries: Sequence[AliasEntry], include_keys_in_file: Set[str], aggressive_unqualified: bool
) -> Tuple[str, bool, Set[str]]:
    changed = False
    out = text
    used_aliases: Set[str] = set()

    # Build a best-match table for base type -> alias to avoid ambiguous
    # mappings (e.g. PatchHierarchy should map to SAMRAIPatchHierarchy, not
    # SAMRAIBasePatchHierarchy).
    best_for_base: Dict[str, Tuple[str, Set[str], bool, int]] = {}
    for entry in entries:
        for base in entry.base_types:
            alias_stem = entry.alias_name[len("SAMRAI") :] if entry.alias_name.startswith("SAMRAI") else entry.alias_name
            score = 3 if alias_stem == base else (2 if alias_stem.endswith(base) else 1)
            prev = best_for_base.get(base)
            if prev is None or score > prev[3]:
                best_for_base[base] = (entry.alias_name, entry.include_keys, entry.is_template_alias, score)

    # Process longer base names first to avoid partial overlaps.
    work: List[Tuple[str, str, Set[str], bool]] = []
    for base, (alias_name, include_keys, is_template_alias, _) in best_for_base.items():
        work.append((base, alias_name, include_keys, is_template_alias))
    work.sort(key=lambda x: len(x[0]), reverse=True)

    for base, alias_name, include_keys, is_template_alias in work:
        # Fully qualified SAMRAI names.
        pat_fq = re.compile(rf"\bSAMRAI::(?:[A-Za-z_][A-Za-z0-9_]*::)+{re.escape(base)}\b")
        out, n1 = pat_fq.subn(alias_name, out)
        if n1:
            changed = True
            used_aliases.add(alias_name)

        # Module-qualified names.
        pat_mod = re.compile(rf"\b(?:{SAMRAI_MODULES_PATTERN})::{re.escape(base)}\b")
        out, n2 = pat_mod.subn(alias_name, out)
        if n2:
            changed = True
            used_aliases.add(alias_name)

        use_unqualified = aggressive_unqualified or base in {"Pointer", "ConstPointer"}
        if not use_unqualified:
            has_qualified_use = _contains_qualified_type(out, base)
            if base in AMBIGUOUS_UNQUALIFIED_BASES and not has_qualified_use:
                continue
            include_hit = any(k in include_keys_in_file for k in include_keys)
            use_unqualified = (
                include_hit
                or has_qualified_use
                or _contains_ndim_template_use(out, base)
            )
        if use_unqualified:
            # Avoid replacing part of a longer identifier (e.g. SAMRAIIntVector).
            pat_unq = re.compile(rf"(?<![A-Za-z0-9_:]){re.escape(base)}\b")
            out, n3 = pat_unq.subn(alias_name, out)
            if n3:
                changed = True
                used_aliases.add(alias_name)

    alias_template_flags = {entry.alias_name: entry.is_template_alias for entry in entries}
    normalized, norm_changed = _normalize_alias_templates(out, alias_template_flags)
    if norm_changed:
        out = normalized
        changed = True

    return out, changed, used_aliases


def _mask_include_lines(text: str) -> Tuple[str, Dict[str, str]]:
    token_to_line: Dict[str, str] = {}
    out_lines: List[str] = []
    for idx, line in enumerate(text.splitlines(keepends=True)):
        line_no_nl = line[:-1] if line.endswith("\n") else line
        if INCLUDE_RE.match(line_no_nl):
            token = f"__SAMRAI_REWRITE_INCLUDE_LINE_{idx}__"
            token_to_line[token] = line
            out_lines.append(token + ("\n" if line.endswith("\n") else ""))
        else:
            out_lines.append(line)
    return "".join(out_lines), token_to_line


def _unmask_include_lines(text: str, token_to_line: Dict[str, str]) -> str:
    out = text
    for token, line in token_to_line.items():
        out = out.replace(token, line[:-1] if line.endswith("\n") else line)
    return out


def _mask_forward_declarations(text: str) -> Tuple[str, Dict[str, str]]:
    token_to_line: Dict[str, str] = {}
    out_lines: List[str] = []
    for idx, line in enumerate(text.splitlines(keepends=True)):
        line_no_nl = line[:-1] if line.endswith("\n") else line
        if FORWARD_DECL_RE.match(line_no_nl):
            token = f"__SAMRAI_REWRITE_FWD_DECL_{idx}__"
            token_to_line[token] = line
            out_lines.append(token + ("\n" if line.endswith("\n") else ""))
        else:
            out_lines.append(line)
    return "".join(out_lines), token_to_line


def _unmask_forward_declarations(text: str, token_to_line: Dict[str, str]) -> str:
    out = text
    for token, line in token_to_line.items():
        out = out.replace(token, line[:-1] if line.endswith("\n") else line)
    return out


def _mask_using_alias_lines(text: str) -> Tuple[str, Dict[str, str]]:
    token_to_line: Dict[str, str] = {}
    out_lines: List[str] = []
    for idx, line in enumerate(text.splitlines(keepends=True)):
        line_no_nl = line[:-1] if line.endswith("\n") else line
        if USING_ALIAS_LINE_RE.match(line_no_nl):
            token = f"__SAMRAI_REWRITE_USING_ALIAS_{idx}__"
            token_to_line[token] = line
            out_lines.append(token + ("\n" if line.endswith("\n") else ""))
        else:
            out_lines.append(line)
    return "".join(out_lines), token_to_line


def _unmask_using_alias_lines(text: str, token_to_line: Dict[str, str]) -> str:
    out = text
    for token, line in token_to_line.items():
        out = out.replace(token, line[:-1] if line.endswith("\n") else line)
    return out


def _mask_template_instantiations(text: str) -> Tuple[str, Dict[str, str]]:
    token_to_line: Dict[str, str] = {}
    out_lines: List[str] = []
    for idx, line in enumerate(text.splitlines(keepends=True)):
        line_no_nl = line[:-1] if line.endswith("\n") else line
        if TEMPLATE_INSTANTIATION_RE.match(line_no_nl):
            token = f"__SAMRAI_REWRITE_TEMPLATE_INST_{idx}__"
            token_to_line[token] = line
            out_lines.append(token + ("\n" if line.endswith("\n") else ""))
        else:
            out_lines.append(line)
    return "".join(out_lines), token_to_line


def _unmask_template_instantiations(text: str, token_to_line: Dict[str, str]) -> str:
    out = text
    for token, line in token_to_line.items():
        out = out.replace(token, line[:-1] if line.endswith("\n") else line)
    return out


def _mask_non_code_regions(text: str) -> Tuple[str, Dict[str, str]]:
    token_to_chunk: Dict[str, str] = {}
    out_parts: List[str] = []
    cursor = 0
    count = 0
    for m in NON_CODE_RE.finditer(text):
        out_parts.append(text[cursor : m.start()])
        token = f"__SAMRAI_REWRITE_NON_CODE_{count}__"
        token_to_chunk[token] = m.group(0)
        out_parts.append(token)
        cursor = m.end()
        count += 1
    out_parts.append(text[cursor:])
    return "".join(out_parts), token_to_chunk


def _unmask_non_code_regions(text: str, token_to_chunk: Dict[str, str]) -> str:
    out = text
    for token, chunk in token_to_chunk.items():
        out = out.replace(token, chunk)
    return out


def _mask_protected_code_tokens(text: str) -> Tuple[str, Dict[str, str]]:
    token_to_chunk: Dict[str, str] = {}
    out = text
    for i, chunk in enumerate(PROTECTED_CODE_TOKENS):
        token = f"__SAMRAI_REWRITE_PROTECTED_{i}__"
        if chunk in out:
            out = out.replace(chunk, token)
            token_to_chunk[token] = chunk
    return out, token_to_chunk


def _unmask_protected_code_tokens(text: str, token_to_chunk: Dict[str, str]) -> str:
    out = text
    for token, chunk in token_to_chunk.items():
        out = out.replace(token, chunk)
    return out


def _split_top_level_template_args(arg_text: str) -> List[str]:
    args: List[str] = []
    depth = 0
    start = 0
    for i, ch in enumerate(arg_text):
        if ch == "<":
            depth += 1
        elif ch == ">":
            depth -= 1
        elif ch == "," and depth == 0:
            args.append(arg_text[start:i].strip())
            start = i + 1
    tail = arg_text[start:].strip()
    if tail:
        args.append(tail)
    return args


def _normalize_alias_templates(text: str, alias_template_flags: Dict[str, bool]) -> Tuple[str, bool]:
    changed = False
    i = 0
    n = len(text)
    out: List[str] = []
    token_re = re.compile(r"\b(SAMRAI[A-Za-z0-9_]+)\b")
    while i < n:
        m = token_re.search(text, i)
        if m is None:
            out.append(text[i:])
            break
        alias_name = m.group(1)
        out.append(text[i : m.end()])
        j = m.end()
        while j < n and text[j].isspace():
            j += 1
        if alias_name in alias_template_flags and j < n and text[j] == "<":
            depth = 0
            k = j
            while k < n:
                if text[k] == "<":
                    depth += 1
                elif text[k] == ">":
                    depth -= 1
                    if depth == 0:
                        k += 1
                        break
                k += 1
            if depth == 0 and k > j:
                template_payload = text[j + 1 : k - 1]
                is_template_alias = alias_template_flags[alias_name]
                if is_template_alias:
                    args = _split_top_level_template_args(template_payload)
                    if args and args[0] == "NDIM":
                        if len(args) > 1:
                            out.append("<" + ", ".join(args[1:]) + ">")
                        changed = True
                        i = k
                        continue
                else:
                    changed = True
                    i = k
                    continue
        i = m.end()
    return "".join(out), changed


def _ensure_names_header(text: str) -> Tuple[str, bool]:
    names_include = '#include "ibtk/samrai_compatibility_names.h"'
    if names_include in text:
        return text, False

    if not re.search(r"\bSAMRAI[A-Z][A-Za-z0-9_]*\b", text):
        return text, False

    lines = text.splitlines()
    insert_at = None
    for i, line in enumerate(lines):
        if line.strip() == "#include <ibtk/config.h>":
            insert_at = i + 1
            break
    if insert_at is None:
        for i, line in enumerate(lines):
            if line.strip().startswith("#include "):
                insert_at = i
                break
    if insert_at is None:
        insert_at = 0
    lines.insert(insert_at, names_include)
    return "\n".join(lines) + ("\n" if text.endswith("\n") else ""), True


def _ensure_alias_includes(text: str, alias_includes: Set[str]) -> Tuple[str, bool]:
    if not alias_includes:
        return text, False

    existing: Set[str] = set()
    lines = text.splitlines()
    for line in lines:
        m = INCLUDE_RE.match(line)
        if m:
            existing.add(m.group(3).strip())

    missing = sorted([inc for inc in alias_includes if inc not in existing])
    if not missing:
        return text, False

    insert_at = None
    for i, line in enumerate(lines):
        if line.strip() == "#include <ibtk/samrai_compatibility_names.h>":
            insert_at = i + 1
            break
    if insert_at is None:
        for i, line in enumerate(lines):
            if line.strip() == "#include <ibtk/config.h>":
                insert_at = i + 1
                break
    if insert_at is None:
        for i, line in enumerate(lines):
            if line.strip().startswith("#include "):
                insert_at = i
                break
    if insert_at is None:
        insert_at = 0

    for inc in missing:
        lines.insert(insert_at, f'#include "{inc}"')
        insert_at += 1

    return "\n".join(lines) + ("\n" if text.endswith("\n") else ""), True


def _collect_alias_includes_from_text(text: str, alias_to_include: Dict[str, str]) -> Set[str]:
    # Only scan code regions so aliases mentioned in comments/strings do not
    # affect include insertion.
    masked_text, _ = _mask_non_code_regions(text)
    includes: Set[str] = set()
    for m in re.finditer(r"\b(?:IBTK::)?(SAMRAI[A-Za-z0-9_]+)\b", masked_text):
        token = m.group(1)
        inc = alias_to_include.get(token)
        if inc is not None:
            includes.add(inc)
    return includes


def _remove_ibtk_alias_qualification(text: str, alias_names: Set[str]) -> Tuple[str, bool]:
    changed = False

    def _repl(match: re.Match[str]) -> str:
        nonlocal changed
        token = match.group(1)
        if token in alias_names:
            changed = True
            return token
        return match.group(0)

    masked_text, non_code_tokens = _mask_non_code_regions(text)
    new_text = re.sub(r"\bIBTK::(SAMRAI[A-Za-z0-9_]+)\b", _repl, masked_text)
    new_text = _unmask_non_code_regions(new_text, non_code_tokens)
    return new_text, changed


def transform_file(
    path: pathlib.Path, repo_root: pathlib.Path, entries: Sequence[AliasEntry], aggressive_unqualified: bool
) -> Tuple[str, str]:
    old_text = path.read_text(encoding="utf-8")
    old_lines = old_text.splitlines()
    new_lines, include_info, includes_changed = _replace_includes(old_lines, entries)
    text_after_includes = "\n".join(new_lines) + ("\n" if old_text.endswith("\n") else "")
    masked_text, include_tokens = _mask_include_lines(text_after_includes)
    masked_text, fwd_decl_tokens = _mask_forward_declarations(masked_text)
    masked_text, using_alias_tokens = _mask_using_alias_lines(masked_text)
    masked_text, template_inst_tokens = _mask_template_instantiations(masked_text)
    masked_text, non_code_tokens = _mask_non_code_regions(masked_text)
    masked_text, protected_tokens = _mask_protected_code_tokens(masked_text)

    text_after_types, types_changed, used_aliases = _replace_types(
        masked_text, entries, include_info["file_include_keys"], aggressive_unqualified
    )

    text_after_types = _unmask_protected_code_tokens(text_after_types, protected_tokens)
    text_after_types = _unmask_non_code_regions(text_after_types, non_code_tokens)
    text_after_types = _unmask_template_instantiations(text_after_types, template_inst_tokens)
    text_after_types = _unmask_using_alias_lines(text_after_types, using_alias_tokens)
    text_after_types = _unmask_forward_declarations(text_after_types, fwd_decl_tokens)
    text_after_types = _unmask_include_lines(text_after_types, include_tokens)
    alias_names = {entry.alias_name for entry in entries}
    text_after_types, dequalify_changed = _remove_ibtk_alias_qualification(text_after_types, alias_names)

    alias_to_include = {entry.alias_name: entry.include_name for entry in entries}
    used_alias_includes = {alias_to_include[a] for a in used_aliases if a in alias_to_include}
    used_alias_includes.update(_collect_alias_includes_from_text(text_after_types, alias_to_include))
    text_with_alias_includes, alias_includes_changed = _ensure_alias_includes(text_after_types, used_alias_includes)
    text_with_names, names_changed = _ensure_names_header(text_with_alias_includes)

    if includes_changed or types_changed or dequalify_changed or alias_includes_changed or names_changed:
        return old_text, text_with_names
    return old_text, old_text


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", nargs="+", help="File(s) to rewrite.")
    parser.add_argument("--write", action="store_true", help="Write changes in place. Default is dry-run diff.")
    parser.add_argument(
        "--aggressive-unqualified",
        action="store_true",
        help="Always rewrite unqualified SAMRAI base names (higher collision risk).",
    )
    args = parser.parse_args(argv)

    repo_root = _repo_root(pathlib.Path(__file__))
    compat_root = (_compat_dir(repo_root)).resolve()
    entries = _load_alias_catalog(repo_root)
    if not entries:
        print("error: no SAMRAI compatibility alias entries found", file=sys.stderr)
        return 2

    changed_files = 0
    for file_arg in args.file:
        path = pathlib.Path(file_arg).resolve()
        try:
            path.relative_to(compat_root)
            print(f"skipped (compat header definition file): {path}")
            continue
        except ValueError:
            pass
        old_text, new_text = transform_file(path, repo_root, entries, args.aggressive_unqualified)
        if old_text == new_text:
            print(f"unchanged: {path}")
            continue
        changed_files += 1
        if args.write:
            path.write_text(new_text, encoding="utf-8")
            print(f"updated: {path}")
        else:
            diff = difflib.unified_diff(
                old_text.splitlines(),
                new_text.splitlines(),
                fromfile=str(path),
                tofile=str(path),
                lineterm="",
            )
            print("\n".join(diff))

    if not args.write and changed_files == 0:
        print("no changes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
