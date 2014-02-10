#!/bin/bash

./writetest
cmp -s writetest.out writetest.checkL || cmp -s writetest.out writetest.checkB
