#!/bin/sh
now=$(date +"%m_%d_%Y")
file="elmerfem-tests-$now.tar.gz"
echo "Packing tests to file: $file..."
tar -czvf $file tests
