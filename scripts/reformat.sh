#! /bin/bash
FILES="`find src -name *.cpp -o -name *.h` `find ibtk/src -name *.cpp -o -name *.h` `find examples -name *.cpp -o -name *.h` `find ibtk/examples -name *.cpp -o -name *.h`"
for file in $FILES; do
  echo $file
  clang-format -i -style=file $file
done
