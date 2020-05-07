#!/usr/bin/env bash
SRC_DIR=$(pwd)
NJOB=${1:-1} # set to 1 if 1st command line argument is not given
KEY=${2:-Release} # set to Release if argument is not given

declare -a EXTRA
if [[ $KEY == "Debug" ]]; then
    BT=Debug
    BUILD_DIR=../cmake_debug
    # EXTRA=(-DCMAKE_CXX_FLAGS="-fsanitize=address -g" -DCMAKE_C_FLAGS="-fsanitize=address -g" -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address -lasan" -DCMAKE_MODULE_LINKER_FLAGS="-fsanitize=address -lasan")
elif [[ $KEY == "Release" ]]; then
    BT=Release
    BUILD_DIR=../cmake_release
else
    echo Build type $KEY not recognized.
    exit 0
fi

mkdir -p $BUILD_DIR
cd $BUILD_DIR

cmake $SRC_DIR -DCMAKE_BUILD_TYPE=$BT \
      "${EXTRA[@]}" \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=YES -DCMAKE_INSTALL_PREFIX=$BUILD_DIR

make -j$NJOB
make install
ln -sf $BUILD_DIR/compile_commands.json $SRC_DIR
