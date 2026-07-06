#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

echo "==========================================================="
echo " Building External Dependencies (htslib + bzip2)"
echo "==========================================================="

mkdir -p external/deps/lib
mkdir -p external/deps/include

# 1. Build bzip2
if [ ! -f "external/deps/lib/libbz2.a" ]; then
    echo "Building bzip2..."
    cd external/deps
    # Since we can't fetch easily inside CI without curl/wget issues, we will just use 
    # the system package managers on macOS/Ubuntu if we were running in CI.
    # But wait, this script compiles it from source for universal compatibility!
    
    # Download bzip2 source
    curl -sLO https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz
    tar -xzf bzip2-1.0.8.tar.gz
    cd bzip2-1.0.8
    make -j4 libbz2.a > /dev/null 2>&1
    cp libbz2.a ../lib/
    cp bzlib.h ../include/
    cd ../..
    rm -rf deps/bzip2-1.0.8*
    echo "bzip2 built successfully."
else
    echo "bzip2 already built."
fi

# 2. Build htslib
if [ ! -f "external/htslib/libhts.a" ]; then
    echo "Building htslib..."
    cd external/htslib
    # htslib requires autoreconf
    autoreconf -i > /dev/null 2>&1
    ./configure --disable-bz2 --disable-lzma --disable-libcurl > /dev/null 2>&1
    make clean > /dev/null 2>&1
    make lib-static NONCONFIGURE_OBJS="" CFLAGS="-g -O2 -fPIC -I$(pwd)/../deps/include" -j4 > /dev/null 2>&1
    cd ../..
    echo "htslib built successfully."
else
    echo "htslib already built."
fi

echo "All dependencies ready!"
