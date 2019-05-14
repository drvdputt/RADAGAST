SRC_DIR=$(pwd)
NJOB=${1:-1} # set to 1 if 1st command line argument is not given
BUILD_DIR=../cmake_release
mkdir $BUILD_DIR
cd $BUILD_DIR

cmake -DCMAKE_BUILD_TYPE=Release $SRC_DIR
ln -s compile_commands.json $SRC_DIR
make -j$NJOB
