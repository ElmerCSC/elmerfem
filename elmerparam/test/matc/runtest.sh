#!/bin/bash

rm -f input

cat > input << STOP
6
1.0
2.0
3.0
4.0
5.0
6.0
6
0
1
2
3
4
5
STOP

ElmerParam input | awk 'NF == 1{if ($1 == 6.0) exit(0); else exit(1)}'
