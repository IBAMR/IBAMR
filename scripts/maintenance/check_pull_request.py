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
        [--title TITLE] [--body FILE]

checks the commits on the current branch that are not on REF (default:
origin/master; for a stacked branch, give the branch it is based on). With
--title and --body it also checks the intended title and description.

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

This is a screen, not a review. It cannot judge whether a comment is clear or
whether a test case is needed, and a reported line may be legitimate. The exit
status is nonzero when anything other than the line counts is reported.
"""

import argparse
import collections
import hashlib
import json
import re
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


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--base", default="origin/master")
    parser.add_argument("--title")
    parser.add_argument("--body", help="file containing the description")
    parser.add_argument("--pr", type=int, help="number of an existing PR")
    parser.add_argument("--repo", default="IBAMR/IBAMR")
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
    failed = False
    for section, notes in (
        ("description", check_description(title, body, labels, base)),
        ("tests", test_notes),
        ("comments", check_comments(files)),
    ):
        for note in notes:
            print(f"{section}: {note}")
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
