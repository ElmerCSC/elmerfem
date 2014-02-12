#!/bin/bash

make clean && make

rm -f save.dat
./foobar > result
