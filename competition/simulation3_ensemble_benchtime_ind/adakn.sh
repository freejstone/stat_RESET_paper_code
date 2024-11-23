#!/bin/bash

values=($(seq 0 4))

for value in "${values[@]}"; do
    Rscript simulation3_ensemble_benchtime_ind/adakn.R $value 10000 3
done