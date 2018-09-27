# Changelog entries.

This folder contains short log entries that describe changes to IBAMR. We will
use this information to generate summaries describing each release of IBAMR: new
features, incompatibilities, etc. At the present time we have the following
categories:

* `major`
* `minor`
* `incompatibilities`

Each log entry should contain the names of the authors, the date, and a short
description of the change, e.g.,
```
New: We can now use the nine-point stencil.
<br>
(John Doe, YYYY/MM/DD)
```
which is named `YYYYMMDDJohnDoe`. File names for multiple contributions from the
same author can have a number appended, e.g., `YYYYMMDDJohnDoe_1`. The file name
will be used to correctly order changelog entries in a generated file at the
time of the next release.

## incompatibilities
This section includes entries that describe incompatible changes to the library
(e.g., requiring C++11 or changing an API).

## major
This section includes major changes, e.g., adding a new tutorial program.

## minor
Finally, the "minor" section contains all the changes that do not fit in
the former two. Bug fixes and pull requests that generalize functions
are typical examples.
