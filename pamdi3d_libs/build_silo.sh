#!/usr/bin/env bash

# Set names/directories
TARGET_DIR=`pwd`"/silo_lib"
SILO_BASEURL="https://wci.llnl.gov/content/assets/docs/simulation/"\
"computer-codes/silo/silo-4.10/"
SILO_TARNAME="silo-4.10-bsd-smalltest.tar.gz"
SILO_DIRNAME="silo-4.10-bsd"
BUILD_DIR="build"

# Do compilation etc. in build directory
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# Get silo if not found
if [ ! -f ${SILO_TARNAME} ]; then
    if hash curl 2> /dev/null; then
        curl -O ${SILO_BASEURL}${SILO_TARNAME}
    elif hash wget 2> /dev/null; then
        wget ${SILO_BASEURL}${SILO_TARNAME}
    else
        echo "build_silo.sh error: cannot find curl or wget"
    fi
fi

# Extract
if [ ! -d ${SILO_DIRNAME} ]; then
    tar -xzf ${SILO_TARNAME}
fi

# Configure
cd ${SILO_DIRNAME}
./configure --disable-shared --disable-fpzip --disable-hzip --disable-silex \
    --disable-browser --disable-dependency-tracking --enable-optimization \
    --disable-libtool-lock --prefix=${TARGET_DIR} --without-hdf5
make
make install
