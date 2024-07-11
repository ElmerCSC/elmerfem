#!/bin/bash

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to install packages
install_packages() {
    echo "Installing required packages..."
    pacman -Syu --noconfirm
    pacman -S --noconfirm mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake git mingw-w64-x86_64-openblas
    pacman -S --noconfirm mingw-w64-x86_64-suitesparse
}

# Check if running in the correct environment
if [[ "$(echo $MSYSTEM)" != "MINGW64" ]]; then
    echo "Error: This script must be run in the MSYS2 MinGW 64-bit environment."
    exit 1
fi

if [[ ! "$PATH" =~ /mingw64/bin ]]; then
    echo "Error: /mingw64/bin is not in the PATH. Make sure you're using the correct MSYS2 environment."
    exit 1
fi

if [[ "$(gcc -dumpmachine)" != "x86_64-w64-mingw32" ]]; then
    echo "Error: Incorrect GCC version. Expected x86_64-w64-mingw32."
    exit 1
fi

# Check and install required tools
for tool in cmake gfortran gcc; do
    if ! command_exists $tool; then
        echo "$tool is not installed. Installing required packages..."
        install_packages
        break
    fi
done

# Check versions
cmake_version=$(cmake --version | head -n1 | awk '{print $3}')
gfortran_version=$(gfortran --version | head -n1 | awk '{print $NF}')
gcc_version=$(gcc --version | head -n1 | awk '{print $NF}')

echo "CMAKE version: $cmake_version"
echo "GFortran version: $gfortran_version"
echo "GCC version: $gcc_version"

if [[ -z "$cmake_version" || -z "$gfortran_version" || -z "$gcc_version" ]]; then
  echo "Error: One or more required tools (CMake, GFortran, GCC) are not properly installed."
  exit 1
fi

# Check if OpenBLAS is installed
if ! pacman -Qs mingw-w64-x86_64-openblas > /dev/null; then
    echo "OpenBLAS is not installed. Installing..."
    pacman -S --noconfirm mingw-w64-x86_64-openblas
fi

# Create or clean build directory
if [[ ! -d "build" ]]; then
    mkdir build
elif [[ "$1" == "clean" ]]; then
    echo "Cleaning build directory..."
    rm -rf build/*
fi

# Build process
cd build || exit 1
cmake -G "MinGW Makefiles" -DWITH_MPI=OFF ..
if [ $? -ne 0 ]; then
    echo "CMake configuration failed."
    exit 1
fi

mingw32-make 2>&1 | tee build_log.txt
if [ $? -ne 0 ]; then
    echo "Build failed. Check build_log.txt for details."
    exit 1
fi

echo "Build process completed successfully."