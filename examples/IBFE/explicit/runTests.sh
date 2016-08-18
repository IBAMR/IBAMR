#!/bin/bash
#set -e 
#set -o pipefail

for dir in ./*/
do
    dir=${dir%*/}
    cd ${dir##*/}
    make gtest
    if [ -f "test2d" ];
    then
        echo "************Running "test2d"************"
        ./test2d input2d.test --gtest_filter=*.2d
    fi
    if [ -f "test3d" ];
    then
        echo "************Running "test3d"************"
        ./test3d input3d.test --gtest_filter=*.3d
    fi
    cd ..
done
