#!/bin/bash

values=($(seq 0 99))

for value in "${values[@]}"; do
    Rscript simulation2_vary_ensemble/RESET.R $value 
done

for value in "${values[@]}"; do
    Rscript simulation2_vary_ensemble/adakn.R $value
done