#!/bin/bash

n=0
m=0

export PATH="../../src/:${PATH}"

for t in *; do
    if [ -d $t ] && [ $t != CVS ]; then \
        cd $t;
        printf "Testing %-20s" ${t}" ..."
        if ./runtest.sh > test.log 2>&1 ; then
            n=$(($n+1))
            printf "\t[PASSED]\n"
        else
            m=$(($m+1))
            printf "\t[FAILED]\n" $t
        fi
        cd ..
    fi
done
echo 'Successful tests:' $n
echo 'Failed tests:' $m

[ $m -gt 0 ] && exit 1 || exit 0;
