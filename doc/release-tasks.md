# IBAMR release tasks

## Summary

This document describes the process for creating a new release of IBAMR on
GitHub.

## static analysis

- [ ] clear basic compilation warnings with GCC and clang on linux
- [ ] compile with the latest libMesh release (which should be configured with
  `--disable-deprecated`) to verify that we do not use any deprecated
  functionality (this will cause compilation errors if we do use them). Add
  version checks: for example, here is one (lightly edited) in
  `IBFEInstrumentPannel.cpp`:
```cpp
    // new API in 1.4.0
#if 1 <= LIBMESH_MAJOR_VERSION && 4 <= LIBMESH_MINOR_VERSION
    const auto node_list = boundary_info.build_node_list();
    for (const auto& pair : node_list)
    {
        nodes.push_back(std::get<0>(pair));
        bcs.push_back(std::get<1>(pair));
    }
#else
    boundary_info.build_node_list(nodes, bcs);
#endif
```
  around APIs to keep backwards compatibility with older versions of libMesh.
- [ ] check the test suite with the oldest (3.7.0) and newest (latest release)
  versions of PETSc. Also check with the oldest (1.1.0) and newest (latest
  release) versions of libMesh.
- [ ] clear basic compilation warnings with GCC and clang on linux and macOS
- [ ] run cppcheck via
```
  cppcheck --force -j4 --enable=all -I./include/ -I./ibtk/include/  \
  $(find ./src ./ibtk/src -name '*.cpp') >cppcheck-results.txt 2>&1
```
- [ ] Check that we don't call functions like
  `ReplicatedMesh::active_local_elements_end()` inside `for`-loop declarations
  since the end iterator is expensive to compute.
- [ ] Run *include what you use* and delete unnecessary header inclusions.
- [ ] Run clang-format on the entire code base.

## testing

- [ ] Send out a message to the mailing list about the release asking about the
  status of any open PRs or bugs.
- [ ] Run the test suite on a cluster where libMesh and PETSc were built without
  HDF5 support and HDF5 was installed as a user (see #687). Make sure all
  examples can compile and link.
- [ ] Check that everything works on macOS where dependencies are installed with
  HomeBrew.
- [ ] In general, run the test suite on as many platforms and configurations as
  possible to weed out problems.

## Copyright

- [ ] Ensure that the top level file `COPYRIGHT` has the release year in it.
- [ ] Check that new files that have been added since the last release have the
  correct copyright header. Try, for example,
```
  git log --diff-filter=A --numstat --pretty='COMMIT: %H %cd' v0.5.1..master
```
  to print all added files between `0.5.1` and `master`.
- [ ] Run a script that updates (based on commits) the second year in the
  copyright statement of each file: i.e., if a file has been changed in the
  current year, then the copyright should be extended to include this year. Make
  sure that headers, source files, automake files, etc. are all updated.

## creating the release

- [ ] Collate the changelog entries on `master`. At the current time this is
  done manually. Open a PR for this on GitHub: only proceed once it is merged.
- [ ] Create a new branch `IBAMR-X.Y` locally off of master and increment the
  version in `VERSION` correctly. Push to GitHub without a PR:
```
  git push https://github.com/ibamr/ibamr.git IBAMR-X.Y
```
  where `X` and `Y` are the version numbers for the release (e.g., for the
  `0.5.0` release they are `0` and `5`). Note that the bugfix number does not go
  into the branch name.
- [ ] Increment the version number on `master` to, e.g., `0.6.0-pre` to signify
  that `master` now corresponds to a prerelease state. Put this in a PR.
- [ ] At this point the release branch and `master` will have diverged. Create a
  release candidate from the `IBAMR-X.Y` branch on GitHub and tag it as a
  prerelease (the tag should be `vX.Y.Z-rc1`). Be sure to include instructions
  on how to patch and compile SAMRAI.
- [ ] If necessary, patch `IBAMR-X.Y` should the release candidate uncover more
  bugs.
- [ ] Once `IBAMR-X.Y` is in good shape tag a proper release, i.e., `vX.Y.Z` and
  create a non-draft release on GitHub.
- [ ] Send out a message to the mailing list about the release.
- [ ] Add a note on the release to the web page.
- [ ] Update version numbers on the build instructions web page.
- [ ] Celebrate!
