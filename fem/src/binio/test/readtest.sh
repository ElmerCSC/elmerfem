#!/bin/bash
ok=`./readtest|sed -e "s/^ *//`
[ "$ok" = " OK" ]