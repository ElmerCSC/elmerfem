#!/bin/bash

cat > input << STOP
3
1.0
2.0
3.0
0
3
STOP

cat input | ElmerParam - output "muu"

if awk '{if ($1 = "muu" && $2 == 1.0 && $3 == 2.0 && $4 == 3.0) exit(1); else exit(0)}' muu; then
    exit 1
fi 

if awk '{if ($1 != NR) exit(0); n++}END{if (n != 3) exit(0); else exit(1)}' output; then
    exit 1
fi

rm -rf input output muu
