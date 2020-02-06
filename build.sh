#!/usr/bin/env bash
SRC_DIR=$(pwd)
NJOB=${1:-1} # set to 1 if 1st command line argument is not given
BT=${2:-Release} # set to Release if argument is not given

if [[ $BT == "Debug" ]]; then
    BUILD_DIR=../cmake_debug
    EXTRA=-DCMAKE_CXX_FLAGS="-fsanitize=address -g" -DCMAKE_C_FLAGS="-fsanitize=address -g" -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address" -DCMAKE_MODULE_LINKER_FLAGS="-fsanitize=address"
elif [[ $BT == "Release" ]]; then
    BUILD_DIR=../cmake_release
else
    BUILD_DIR=../cmake_$BT
fi

if [[ ! -d $BUILD_DIR ]]; then
    mkdir $BUILD_DIR
fi
cd $BUILD_DIR

# Explicitly unset EIGEN_DIR here, to test find_path
cmake $SRC_DIR \
      $EXTRA \
      -DCMAKE_BUILD_TYPE=$BT -DCMAKE_EXPORT_COMPILE_COMMANDS=YES -DCMAKE_INSTALL_PREFIX=$BUILD_DIR -UEIGEN_DIR

ln -sf $BUILD_DIR/compile_commands.json $SRC_DIR
make -j$NJOB

make install
