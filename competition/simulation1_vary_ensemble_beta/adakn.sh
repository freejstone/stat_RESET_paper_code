#!/bin/bash

values=($(seq 0 19))
ks=(150 300 450)

for k in "${ks[@]}"; do
    for value in "${values[@]}"; do
         Rscript simulation1_vary_ensemble_beta/adakn.R $value $k;
    done
done
