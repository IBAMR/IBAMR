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

"""Screen a branch or pull request against the conventions in doc/agents.

Before publishing, from the top-level directory of IBAMR:

    scripts/maintenance/check_pull_request.py [--base REF]
        [--title TITLE] [--body FILE] [--compiler-check BUILD_DIRECTORY]

checks the commits on the current branch that are not on REF (default:
origin/master; for a stacked branch, give the branch it is based on). With
--title and --body it also checks the intended title and description. With
--compiler-check it also compiles the changed C++ sources, syntax only, with
the compiler family that BUILD_DIRECTORY does not use: GCC if it was configured
with Clang, Clang if it was configured with GCC. The include paths and
definitions come from BUILD_DIRECTORY/compile_commands.json (configure with
-DCMAKE_EXPORT_COMPILE_COMMANDS=ON). The hosted builds use GCC on Linux and
Clang on macOS and treat warnings as errors, and each compiler reports
problems that the other does not. The check is skipped, with a note, when no
such compiler or no compilation database is found.

For a pull request that already exists:

    scripts/maintenance/check_pull_request.py --pr NUMBER

reads the title, description, labels, and diff with the GitHub CLI.

The script reports:

  description  length, test lists, defect history, template, labels, base
  tests        added lines of library code, test code, input, and expected
               output, and expected outputs identical to another added by the
               same change
  comments     added comment lines containing development history or working
               shorthand, and added uses of typeid
  compiler     with --compiler-check, errors and warnings that the other
               compiler reports for the changed C++ sources

This is a screen, not a review. It cannot judge whether a comment is clear or
whether a test case is needed, and a reported line may be legitimate. The exit
status is nonzero when anything other than the line counts is reported.
"""

import argparse
import collections
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys

TEMPLATE_COMMENT = "This template should be included in all pull requests."
TEMPLATE_HEADING = "### IBAMR Pull Request Checklist"
MAX_DESCRIPTION_WORDS = 60
MAX_TITLE_LENGTH = 72

HISTORY = re.compile(
    r"previously|no longer|legacy|retained|used to\b|this fix|"
    r"next PR|in this stack|master's|before this change|\bnow\b",
    re.IGNORECASE,
)
SHORTHAND = re.compile(
    r"\bparity\b|\bpins?\b|\blive\b|\bprobe\b", re.IGNORECASE
)
ATTRIBUTION = re.compile(r"Co-Authored-By|Generated with", re.IGNORECASE)
FIX_TITLE = re.compile(r"(\[[^\]]*\]\s*)*Fix\b")
DYNAMIC_TYPE = re.compile(r"typeid\s*\(")

SOURCE_SUFFIXES = (".h", ".cpp", ".m4", ".f", ".i")
FORTRAN_SUFFIXES = (".m4", ".f", ".i")


def run(*command):
    return subprocess.run(
        command, check=True, capture_output=True, text=True
    ).stdout


def is_comment(path, text):
    stripped = text.strip()
    if path.endswith(FORTRAN_SUFFIXES):
        # Fixed-form comment lines have c or C in the first column.
        return text[:1] in ("c", "C") or "!" in stripped
    return "//" in stripped or stripped.startswith(("*", "/*"))


def code_part(path, text):
    """Return the part of a C++ line that is not a comment."""
    stripped = text.strip()
    if path.endswith(FORTRAN_SUFFIXES) or stripped.startswith(("*", "/*")):
        return ""
    return stripped.split("//")[0]


def parse_diff(diff):
    """Return {path: {"new": bool, "added": [(line_number, text)]}}."""
    files = collections.OrderedDict()
    path = None
    line_number = 0
    for line in diff.splitlines():
        if line.startswith("diff --git"):
            path = re.search(r" b/(.*)$", line).group(1)
            files[path] = {"new": False, "added": []}
        elif path is None:
            continue
        elif line.startswith("new file mode"):
            files[path]["new"] = True
        elif line.startswith("@@"):
            line_number = int(re.search(r"\+(\d+)", line).group(1))
        elif line.startswith("+") and not line.startswith("+++"):
            files[path]["added"].append((line_number, line[1:]))
            line_number += 1
        elif not line.startswith("-"):
            line_number += 1
    return files


def check_description(title, body, labels, base):
    notes = []
    if title is not None:
        if len(title) > MAX_TITLE_LENGTH:
            notes.append(f"the title is {len(title)} characters long")
    if body is None:
        return notes
    body = body.replace("\r", "")
    opening = body.split("<!--")[0].split(TEMPLATE_HEADING)[0].strip()
    words = len(opening.split())
    if words > MAX_DESCRIPTION_WORDS:
        notes.append(
            f"{words} words precede the template; begin with one or two "
            "sentences and add only what a reviewer needs"
        )
    if re.search(r"^Tests?:", opening, re.M):
        notes.append("the description lists the tests")
    if re.search(r"^Validation:", opening, re.M):
        notes.append("the description contains a validation report")
    if re.search(r"\b[Pp]reviously\b|On current master", opening):
        notes.append("the description narrates the history of the defect")
    if re.match(
        r"(Adds|Fixes|Removes|Corrects|Replaces|Supports)\b", opening
    ) and not re.match(r"Fixes #\d", opening):
        notes.append("the opening sentence is not imperative")
    if ATTRIBUTION.search(body):
        notes.append("the description contains an attribution line")
    if TEMPLATE_HEADING in body and TEMPLATE_COMMENT not in body:
        notes.append("the template is missing its leading comment")
    if base not in (None, "master", "origin/master") and not re.search(
        r"Stacked on #\d+", opening
    ):
        notes.append(f"the base is {base}, but no 'Stacked on #N.' is given")
    if labels is not None:
        if not labels:
            notes.append("no labels are set")
        elif title is not None and ("Bug" in labels) != bool(
            FIX_TITLE.match(title)
        ):
            notes.append("the title and the Bug label disagree")
    return notes


def count_tests(files):
    counts = collections.Counter()
    outputs = collections.defaultdict(list)
    for path, info in files.items():
        added = len(info["added"])
        if path.startswith("tests/"):
            if path.endswith(".input"):
                counts["input"] += added
                counts["new_inputs"] += info["new"]
            elif path.endswith(".output"):
                counts["output"] += added
                # An assertion-only case can have empty expected output.
                if info["new"] and info["added"]:
                    text = "\n".join(text for _, text in info["added"])
                    digest = hashlib.sha256(text.encode()).hexdigest()
                    outputs[digest].append(path)
            elif path.endswith((".cpp", ".h")):
                counts["test"] += added
                counts["new_sources"] += info["new"] and path.endswith(".cpp")
        elif path.endswith(SOURCE_SUFFIXES) and not path.startswith("examples/"):
            counts["library"] += added
    summary = (
        f"added lines: library {counts['library']}, test source "
        f"{counts['test']}, input {counts['input']}, expected output "
        f"{counts['output']}; new input files: {counts['new_inputs']}; new "
        f"test sources: {counts['new_sources']}"
    )
    notes = []
    for paths in outputs.values():
        if len(paths) > 1:
            notes.append(
                f"{len(paths)} identical expected outputs: "
                + ", ".join(sorted(paths))
            )
    return summary, notes


def check_comments(files):
    notes = []
    for path, info in files.items():
        if path.startswith("doc/") or not path.endswith(
            SOURCE_SUFFIXES + (".input",)
        ):
            continue
        in_library = not path.startswith(("tests/", "examples/"))
        for line_number, text in info["added"]:
            excerpt = text.strip()[:80]
            if is_comment(path, text):
                for label, pattern in (
                    ("history", HISTORY),
                    ("shorthand", SHORTHAND),
                ):
                    match = pattern.search(text)
                    if match:
                        notes.append(
                            f'{path}:{line_number}: {label} "{match.group(0)}"'
                            f": {excerpt}"
                        )
            if in_library and DYNAMIC_TYPE.search(code_part(path, text)):
                notes.append(f"{path}:{line_number}: typeid: {excerpt}")
    return notes


# The warning flags of the hosted builds for each compiler family. Warnings are
# reported, not promoted to errors, so that every one is listed.
COMPILER_FLAGS = {
    "GCC": ("-Wall", "-Wextra", "-Wpedantic", "-Wno-deprecated-declarations"),
    "Clang": ("-Wall", "-Wextra", "-Wpedantic", "-Wmost", "-Wmove", "-Wunused"),
}
COMPILER_NAMES = {
    "GCC": ("g++-15", "g++-14", "g++-13", "g++-12", "g++-11", "g++"),
    "Clang": ("clang++", "clang++-20", "clang++-19", "clang++-18", "clang++-17"),
}


def compiler_family(compiler):
    """Return "GCC", "Clang", or None for a compiler or compiler wrapper.

    The name does not decide: g++ is Clang on macOS, and an MPI wrapper can be
    either.
    """
    path = shutil.which(compiler)
    if path is None:
        return None
    version = subprocess.run(
        [path, "--version"], capture_output=True, text=True
    ).stdout
    if "Free Software Foundation" in version:
        return "GCC"
    if "clang" in version.lower():
        return "Clang"
    return None


def command_tokens(entry):
    if "arguments" in entry:
        return list(entry["arguments"])
    return shlex.split(entry["command"])


def other_compiler(entries):
    """Return (compiler, family) of the family the build does not use."""
    build_family = None
    for token in command_tokens(entries[0])[:3]:
        if token.startswith("-"):
            break
        build_family = compiler_family(token) or build_family
    if build_family is None:
        return None, "the compiler of the build directory is neither GCC nor Clang"
    family = "Clang" if build_family == "GCC" else "GCC"
    for name in COMPILER_NAMES[family]:
        if compiler_family(name) == family:
            return shutil.which(name), family
    return None, f"the build uses {build_family} and no {family} compiler was found"


def compiler_arguments(entry, top, build_directory):
    """Return the definitions and include paths of one compilation command.

    Warning and code generation options are dropped because they were written
    for another compiler. Include directories outside the IBAMR sources and the
    build directory, and bundled third-party sources, become system include
    directories so that only IBAMR code is diagnosed.
    """
    tokens = command_tokens(entry)
    arguments = []
    label = ""
    index = 1
    while index < len(tokens):
        token = tokens[index]
        index += 1
        if token in ("-I", "-isystem", "-D", "-U", "-include"):
            value = tokens[index]
            index += 1
        elif token[:2] in ("-I", "-D", "-U"):
            token, value = token[:2], token[2:]
        elif token.startswith("-std="):
            arguments.append(token)
            continue
        else:
            continue
        if token == "-I":
            path = os.path.normpath(os.path.join(entry["directory"], value))
            ours = (
                path.startswith(top + os.sep) and "contrib" not in path
            ) or path.startswith(build_directory + os.sep)
            arguments += ["-I" if ours else "-isystem", path]
        else:
            arguments += [token, value]
            if token == "-D" and value.startswith("NDIM="):
                label = value
    return arguments, label


def check_compiler(files, build_directory, compiler):
    """Return (notes, reason the check was skipped or None)."""
    build_directory = os.path.abspath(build_directory)
    database = os.path.join(build_directory, "compile_commands.json")
    if not os.path.exists(database):
        return [], f"{database} does not exist"
    top = run("git", "rev-parse", "--show-toplevel").strip()
    sources = {
        os.path.join(top, path): path
        for path in files
        if path.endswith(".cpp") and os.path.exists(os.path.join(top, path))
    }
    if not sources:
        return [], None
    with open(database) as database_file:
        entries = [
            entry for entry in json.load(database_file)
            if os.path.normpath(os.path.join(entry["directory"], entry["file"]))
            in sources
        ]
    if not entries:
        return [], f"{database} lists none of the changed C++ sources"
    if compiler is None:
        compiler, family = other_compiler(entries)
        if compiler is None:
            return [], family
    else:
        family = compiler_family(compiler)
        if family is None:
            return [], f"{compiler} is neither GCC nor Clang"
        compiler = shutil.which(compiler)
    notes = []
    seen = set()
    for entry in entries:
        source = os.path.normpath(os.path.join(entry["directory"], entry["file"]))
        arguments, label = compiler_arguments(entry, top, build_directory)
        if (source, label) in seen:
            continue
        seen.add((source, label))
        result = subprocess.run(
            [compiler, "-fsyntax-only", *COMPILER_FLAGS[family], *arguments, source],
            capture_output=True, text=True,
        )
        for line in result.stderr.splitlines():
            if " error: " in line or " warning: " in line:
                where = f" ({family}, {label})" if label else f" ({family})"
                notes.append(line.replace(top + os.sep, "") + where)
    missing = sorted(set(sources.values()) - {sources[s] for s, _ in seen})
    for path in missing:
        notes.append(f"{path}: not in {database}, so not compiled")
    return notes, None


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--base", default="origin/master")
    parser.add_argument("--title")
    parser.add_argument("--body", help="file containing the description")
    parser.add_argument("--pr", type=int, help="number of an existing PR")
    parser.add_argument("--repo", default="IBAMR/IBAMR")
    parser.add_argument(
        "--compiler-check", metavar="BUILD_DIRECTORY",
        help="also compile the changed C++ sources, syntax only, with the "
        "compiler family the build directory does not use",
    )
    parser.add_argument(
        "--compiler", help="compiler for --compiler-check, in place of the "
        "one found automatically",
    )
    args = parser.parse_args()

    labels = None
    if args.pr is not None:
        pr = json.loads(
            run(
                "gh", "pr", "view", str(args.pr), "--repo", args.repo,
                "--json", "title,body,labels,baseRefName",
            )
        )
        title = pr["title"]
        body = pr["body"] or ""
        labels = [label["name"] for label in pr["labels"]]
        base = pr["baseRefName"]
        diff = run("gh", "pr", "diff", str(args.pr), "--repo", args.repo)
    else:
        title = args.title
        body = None
        if args.body is not None:
            with open(args.body) as body_file:
                body = body_file.read()
        base = args.base
        diff = run("git", "diff", f"{base}...HEAD")

    files = parse_diff(diff)
    summary, test_notes = count_tests(files)
    print(summary)
    compiler_notes = []
    if args.compiler_check is not None:
        if args.pr is not None:
            skipped = "it needs the branch checked out; do not use --pr"
        else:
            compiler_notes, skipped = check_compiler(
                files, args.compiler_check, args.compiler
            )
        if skipped is not None:
            print(f"compiler: skipped: {skipped}")
    failed = False
    for section, notes in (
        ("description", check_description(title, body, labels, base)),
        ("tests", test_notes),
        ("comments", check_comments(files)),
        ("compiler", compiler_notes),
    ):
        for note in notes:
            print(f"{section}: {note}")
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
