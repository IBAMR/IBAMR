#! /bin/bash
FILES="`find examples include ibtk/include ibtk/examples ibtk/src src -name *.cpp -o -name *.h`"
for file in $FILES; do
  echo $file
  clang-format -i -style=file $file
done
