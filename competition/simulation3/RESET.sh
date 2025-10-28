#!/bin/bash

values=($(seq 0 99))

for value in "${values[@]}"; do
    Rscript simulation3/RESET.R $value 
done
