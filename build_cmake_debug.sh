SRC_DIR=$(pwd)
NJOB=${1:-1} # set to 1 if 1st command line argument is not given
BUILD_DIR=../cmake_debug
mkdir $BUILD_DIR
cd $BUILD_DIR

cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=YES $SRC_DIR
ln -s compile_commands.json $SRC_DIR
make -j$NJOB
