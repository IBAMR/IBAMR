#! /bin/bash
FILES="`find src -name *.cpp` `find include/ibamr -name *.h` `find ibtk/src -name *.cpp` `find ibtk/include/ibtk -name *.h`"
for file in $FILES; do
  echo $file
  clang-format -i -style=file $file
done
