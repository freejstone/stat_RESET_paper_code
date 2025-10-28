#!/bin/bash

values=($(seq 0 99))
ks=(50 150 300)

for k in "${ks[@]}"; do
    for value in "${values[@]}"; do
        Rscript simulation1/RESET.R $value $k;
    done
done
