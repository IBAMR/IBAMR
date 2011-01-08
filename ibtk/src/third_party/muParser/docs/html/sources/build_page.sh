#!/bin/bash

#
# This script is used for building the muParser html page from
# html templates.
#

rm -rf ../*.html

#
# add navigation bar to all html templates starting with mup_*
#
for file in mup_*
do
  echo processing $file
  cat navigation.html | sed "/\$PLACEHOLDER/r $file" > ../$file
done

# create index.html
cp ../mup_intro.html ../index.html
cat ../mup_intro.html | sed "/\$STAT_COUNTER/r stat_counter.html" > ../index.html

