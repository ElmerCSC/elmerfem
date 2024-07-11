#!/bin/bash

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to install packages
check_and_install_packages() {
    local packages=("mingw-w64-x86_64-toolchain" "mingw-w64-x86_64-cmake" "git" "mingw-w64-x86_64-openblas" "mingw-w64-x86_64-suitesparse")
    local to_install=()
    local to_update=()

    for package in "${packages[@]}"; do
        if ! pacman -Qi "$package" &> /dev/null; then
            to_install+=("$package")
        elif pacman -Qu "$package" &> /dev/null; then
            to_update+=("$package")
        fi
    done

    if [ ${#to_install[@]} -ne 0 ]; then
        echo "The following packages will be installed:"
        printf '%s\n' "${to_install[@]}"
        pacman -S --noconfirm "${to_install[@]}"
    fi

    if [ ${#to_update[@]} -ne 0 ]; then
        echo "The following packages will be updated:"
        printf '%s\n' "${to_update[@]}"
        pacman -S --noconfirm "${to_update[@]}"
    fi

    if [ ${#to_install[@]} -eq 0 ] && [ ${#to_update[@]} -eq 0 ]; then
        echo "All required packages are already installed and up to date."
    fi
}

# Check if running in the correct environment
if [[ "$(echo $MSYSTEM)" != "MINGW64" ]]; then
    echo "Error: This script must be run in the MSYS2 MinGW 64-bit environment."
    exit 1
fi

# Check and install required packages
# check_and_install_packages

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

# Run CMake without UMFPACK_ROOT for now
cmake -G "MinGW Makefiles" -DWITH_MPI=OFF ..
cmake_exit_code=$?

if [ $cmake_exit_code -ne 0 ]; then
    echo "CMake configuration failed."
    exit 1
fi

# Check if UMFPACK was found
if ! grep -q "UMFPACK found:" CMakeCache.txt; then
    echo "Error: UMFPACK was not found by CMake. Attempting to specify UMFPACK_ROOT..."
    
    # Try again with UMFPACK_ROOT specified
    cmake -G "MinGW Makefiles" -DWITH_MPI=OFF -DUMFPACK_ROOT=/mingw64 ..
    cmake_exit_code=$?

    if [ $cmake_exit_code -ne 0 ]; then
        echo "CMake configuration failed even with UMFPACK_ROOT specified."
        exit 1
    fi

    if ! grep -q "UMFPACK found:" CMakeCache.txt; then
        echo "Error: UMFPACK was still not found by CMake. Please check your installation."
        exit 1
    fi
fi

echo "Starting build process..."
mingw32-make 2>&1 | tee build_log.txt
make_exit_code=${PIPESTATUS[0]}

if [ $make_exit_code -ne 0 ]; then
    echo "Build failed. Check build_log.txt for details."
    exit 1
fi

echo "Build process completed successfully."