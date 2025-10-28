#!/bin/bash

values=($(seq 0 4))

for value in "${values[@]}"; do
    Rscript simulation4/adakn.R $value 10000 3
done