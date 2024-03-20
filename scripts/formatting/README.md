# Reformatting Scripts

These scripts reformat source code with a specific clang-format. Afterwards they
perform a few other transformations to ensure a uniform style (e.g., converting
all files to use Unix line endings).

There are three ways to use these scripts:
1. Run `make indent` in the build directory
2. Run `./scripts/format/indent` in the root source directory
3. Run `setup_clang_format_precommit_hook.sh` to set up `indent` to run as a
   pre-commit hook.

# Legal Notice

These scripts are based on deal.II and are licensed under the LGPL - hence they
are not installed with IBAMR nor do the IBAMR developers claim any copyright on
them.
