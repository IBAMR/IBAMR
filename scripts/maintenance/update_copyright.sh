#!/bin/bash
## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2021 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

# Originally based on the deal.II script of the same name with substantial
# modifications made to make it apply to IBAMR and IBTK's needs.

if test ! -d ibtk -o ! -d include -o ! -d examples ; then
  echo "*** This script must be run from the top-level directory of IBAMR."
  exit 1
fi

FILES="
  $(echo attest configure.ac ./ibtk/configure.ac CMakeLists.txt)
  $(find -L ./include ./ibtk/include ./src ./ibtk/src ./cmake ./config ./.github | egrep '\.(m|py|pl|h|am|ac|m4|cpp|i|txt)$')
  $(find -L ./tests ./examples ./ibtk/examples ./scripts | egrep '\.(m|py|pl|h|am|ac|m4|cpp|i)$')
  $(find -L ./doc ./scripts ./m4 ./ibtk/m4 ./lib ./ibtk/lib | egrep '\.(m|py|pl|h|am|ac|m4|cpp|i|txt)$')
"
FILES=$(echo $FILES | xargs realpath | sort -u)

for FILE in $FILES ; do
  # get the last year this file was modified from the git log.
  # we don't want to see patches that just updated the copyright
  # year, so output the dates and log messages of the last 3
  # commits, throw away all that mention both the words
  # "update" and "copyright", and take the year of the first
  # message that remains
  #
  # (it should be enough to look at the last 2 messages since
  # ideally no two successive commits should have updated the
  # copyright year. let's err on the safe side and take the last
  # 3 commits.)
  last_year=`git log -n 3 --date=short --format="format:%cd %s" $FILE | \
             egrep -i -v "update.*copyright|copyright.*update" | \
             head -n 1 | \
             perl -p -e 's/^(\d\d\d\d)-.*/\1/g;'`

  # get the first year this file was modified from the actual
  # file. this may predate the git log if the file was copied
  # from elsewhere
  first_year=`cat $FILE | egrep 'Copyright \(C\) [0-9]{4}' | \
              perl -p -e "s/.*Copyright \(C\) (\d{4}).*/\1/g;"`

  echo "Processing $FILE: ${first_year} - ${last_year}"
  if test ! "${first_year}" = "${last_year}" ; then
    perl -pi -e "s/(Copyright \(c\) \d{4})( - \d{4})?(, \d{4}( - \d{4})?)*.*by the IBAMR developers/\1 - ${last_year} by the IBAMR developers/g;" $FILE
  fi
done
