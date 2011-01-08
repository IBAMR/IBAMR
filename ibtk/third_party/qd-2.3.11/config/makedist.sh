#!/bin/bash
PACKAGE_NAME=qd
MAJOR_VERSION=2
MINOR_VERSION=3
PATCH_LEVEL=$1
if [ -z "$PATCH_LEVEL" ]; then
  echo "Usage: makedist.sh patch-level"
  exit
fi

if [ "$PATCH_LEVEL" = "git" ]; then
  PATCH_LEVEL=git-$(git log --max-count=1 --pretty=oneline | cut -c -8)
fi
VERSION=$MAJOR_VERSION.$MINOR_VERSION.$PATCH_LEVEL

echo "Creating $PACKAGE_NAME-$VERSION distribution..."

DIR=/var/tmp/$PACKAGE_NAME-$$
ORIG_DIR=`pwd`

export CXX=g++
export FC=gfortran

mkdir -p $DIR &&
cp -pr . $DIR &&
cd $DIR &&
mv configure.ac configure.old &&
sed "/^define(\[QD_PATCH_VERSION\]/s/devel/$PATCH_LEVEL/" configure.old >configure.ac && 
rm -f configure.old &&
config/autogen.sh &&
./configure &&
git log --no-merges >ChangeLog &&
make doc &&
make distcheck &&
cp $PACKAGE_NAME-$VERSION.tar.gz $ORIG_DIR &&
rm -rf $DIR

