#!/bin/bash

make clean && make

./runopt > /dev/null 2>&1 && diff evals.dat evals.dat.check
