#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List


ALIAS_SUMMARY_RE = re.compile(r"summary:\s+files=\d+\s+updated=(\d+)\s+unchanged=\d+\s+skipped=\d+\s+failed=\d+")
INCLUDE_CHANGED_RE = re.compile(r"files changed:\s+(\d+)")

DEFAULT_EXTENSIONS = (".h", ".hh", ".hpp", ".hxx", ".c", ".cc", ".cpp", ".cxx", ".inl", ".tcc")


def run_command(cmd: List[str], cwd: Path) -> str:
    print("+ " + " ".join(shlex.quote(c) for c in cmd), flush=True)
    result = subprocess.run(cmd, cwd=str(cwd), check=False, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"command failed with exit code {result.returncode}: {' '.join(cmd)}")
    return (result.stdout or "") + (result.stderr or "")


def chunked(items: List[Path], size: int) -> Iterable[List[Path]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


def parse_alias_updated(output: str) -> int:
    updated = 0
    for line in output.splitlines():
        m = ALIAS_SUMMARY_RE.search(line)
        if m:
            updated += int(m.group(1))
    if "summary:" in output and updated == 0:
        return 0
    if "summary:" not in output:
        raise RuntimeError("could not parse alias rewrite summary output")
    return updated


def parse_include_changed(output: str) -> int:
    changed = 0
    found = False
    for line in output.splitlines():
        m = INCLUDE_CHANGED_RE.search(line)
        if m:
            changed += int(m.group(1))
            found = True
    if not found:
        raise RuntimeError("could not parse include rewrite output (files changed)")
    return changed


def run_alias_pass(repo_root: Path, python: str, alias_script: Path, files: List[Path], write: bool, chunk_size: int) -> int:
    total = 0
    for chunk in chunked(files, chunk_size):
        cmd = [python, str(alias_script)]
        if write:
            cmd.append("--write")
        cmd.append("--summary-only")
        cmd.extend(str(p) for p in chunk)
        out = run_command(cmd, cwd=repo_root)
        total += parse_alias_updated(out)
    return total


def run_include_pass(
    repo_root: Path, python: str, include_script: Path, files: List[Path], write: bool, chunk_size: int
) -> int:
    total = 0
    for chunk in chunked(files, chunk_size):
        cmd = [python, str(include_script)]
        if write:
            cmd.append("--write")
        cmd.extend(str(p) for p in chunk)
        out = run_command(cmd, cwd=repo_root)
        total += parse_include_changed(out)
    return total


def collect_files_from_roots(repo_root: Path, roots: List[str], extensions: tuple[str, ...]) -> List[Path]:
    files: List[Path] = []
    for root in roots:
        abs_root = (repo_root / root).resolve()
        if not abs_root.exists():
            raise RuntimeError(f"root does not exist: {root}")
        for p in abs_root.rglob("*"):
            if p.is_file() and p.suffix in extensions:
                files.append(p.relative_to(repo_root))
    files = sorted(set(files))
    return files


def load_files_from_list(repo_root: Path, list_file: Path) -> List[Path]:
    files = []
    repo_root = repo_root.resolve()
    for line in list_file.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        p = Path(line)
        if p.is_absolute():
            p = p.resolve().relative_to(repo_root)
        files.append(p)
    return sorted(set(files))


def main() -> int:
    parser = argparse.ArgumentParser(description="Run SAMRAI compatibility adoption workflow with convergence checks.")
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--file-list", type=Path, help="Path to newline-delimited target file list.")
    parser.add_argument(
        "--roots",
        nargs="+",
        help="Relative roots to scan when --file-list is not provided.",
    )
    parser.add_argument("--extensions", nargs="+", default=list(DEFAULT_EXTENSIONS))
    parser.add_argument("--python", default="python3")
    parser.add_argument(
        "--alias-script",
        type=Path,
        default=Path("scripts/maintenance/rewrite_samrai_compatibility_file.py"),
    )
    parser.add_argument(
        "--include-script",
        type=Path,
        default=Path("scripts/maintenance/rewrite_include_quotes_to_angle_brackets.py"),
    )
    parser.add_argument("--build-dir", type=Path, default=Path("/Users/boyceg/code/ibamr-objs-dbg"))
    parser.add_argument("--indent-target", default="indent-all")
    parser.add_argument("--max-iterations", type=int, default=10)
    parser.add_argument("--chunk-size", type=int, default=120)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    alias_script = (repo_root / args.alias_script).resolve()
    include_script = (repo_root / args.include_script).resolve()

    if args.file_list:
        files = load_files_from_list(repo_root, args.file_list.resolve())
    else:
        if not args.roots:
            raise RuntimeError("provide --file-list or --roots")
        files = collect_files_from_roots(repo_root, args.roots, tuple(args.extensions))

    if not files:
        print("No files matched; nothing to do.")
        return 0

    print(f"Target files: {len(files)}")

    for i in range(1, args.max_iterations + 1):
        updated = run_alias_pass(repo_root, args.python, alias_script, files, write=True, chunk_size=args.chunk_size)
        print(f"[alias write] iteration {i}: updated={updated}")
        if updated == 0:
            break
    else:
        raise RuntimeError("alias rewrite did not converge")

    for i in range(1, args.max_iterations + 1):
        changed = run_include_pass(repo_root, args.python, include_script, files, write=True, chunk_size=args.chunk_size)
        print(f"[include write] iteration {i}: changed={changed}")
        if changed == 0:
            break
    else:
        raise RuntimeError("include rewrite did not converge")

    # Required gate: rewrites must be stable before indent-all.
    alias_dry = run_alias_pass(repo_root, args.python, alias_script, files, write=False, chunk_size=args.chunk_size)
    include_dry = run_include_pass(repo_root, args.python, include_script, files, write=False, chunk_size=args.chunk_size)
    print(f"[pre-indent verify] alias dry-run updated={alias_dry}")
    print(f"[pre-indent verify] include dry-run changed={include_dry}")

    if alias_dry != 0 or include_dry != 0:
        raise RuntimeError("pre-indent dry-run verification failed")

    run_command(["make", "-C", str(args.build_dir), args.indent_target], cwd=repo_root)

    # Informational only: formatters may re-order include spacing/grouping.
    post_alias_dry = run_alias_pass(repo_root, args.python, alias_script, files, write=False, chunk_size=args.chunk_size)
    post_include_dry = run_include_pass(
        repo_root, args.python, include_script, files, write=False, chunk_size=args.chunk_size
    )
    print(f"[post-indent info] alias dry-run updated={post_alias_dry}")
    print(f"[post-indent info] include dry-run changed={post_include_dry}")

    print("Adoption workflow converged before indent-all.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
