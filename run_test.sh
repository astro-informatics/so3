#!/bin/sh

./bin/c/so3_test 128 4 > timing_test_L0128_N0004.txt
./bin/c/so3_test 256 4 > timing_test_L0256_N0004.txt
./bin/c/so3_test 512 4 > timing_test_L0512_N0004.txt

./bin/c/so3_test 16 16 > timing_test_L0016_N0016.txt
./bin/c/so3_test 32 32 > timing_test_L0032_N0032.txt
./bin/c/so3_test 64 64 > timing_test_L0064_N0064.txt
./bin/c/so3_test 128 128 > timing_test_L0128_N0128.txt


