#!/bin/bash
# Temporary launcher: attest otherwise discards successful subprocess stdout/stderr.
# The real executable receives the same input/options; return its exact exit status.
"$0.real" "$@" > executable-stdout.txt 2> executable-stderr.txt
task_status=$?
printf '%s\n' "$task_status" > executable-exit-status.txt
cat executable-stdout.txt
cat executable-stderr.txt >&2
exit "$task_status"
