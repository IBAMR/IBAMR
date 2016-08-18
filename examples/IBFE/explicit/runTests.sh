#!/bin/bash
#set -e 
#set -o pipefail

testDirs=( ex0 ex1 ex3 ex5 ex6 ex7 ex8 )

for i in "${testDirs[@]}"
do
    dir=$i
    cd ${dir##*/}
    make gtest
    if [ -f "test2d" ];
    then
        echo "************Running "test2d" in $i ************"
        ./test2d input2d.test --gtest_filter=*.2d
    fi
    if [ -f "test3d" ];
    then
        echo "************Running "test3d" in $i ************"
        ./test3d input3d.test --gtest_filter=*.3d
    fi
    cd ..
done

