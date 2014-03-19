#! /bin/bash
FILES="`find src -name *.cpp` `find src -name *.h` `find ibtk/src -name *.cpp` `find ibtk/src -name *.h`"
for file in $FILES; do
  echo $file
  $CLANG_FORMAT -i -style=file $file
done
