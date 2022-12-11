#!/bin/bash
echo "[1/4] Compiling pi calculator..."
g++ -fopenmp pi.cpp -o pi

echo "[2/4] chmod +x pi ..."
chmod +x pi

echo "[3/4] Calculating pi via numerical integration (1-8 threads) ..."
./pi 8 num_int 100000000 num_int_results.csv

echo "[4/4] Calculating pi via Monte-Carlo (1-8 threads) ..."
./pi 8 mc 10000000 mc_results.csv
