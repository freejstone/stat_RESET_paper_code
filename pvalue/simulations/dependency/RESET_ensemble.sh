#!/bin/bash

section1() {
    values=($(seq 1 100))
    for value in "${values[@]}"; do 
        Rscript RESET_ensemble.R data sim1_"${value}".csv result_data/reset_all_sim1_"${value}".csv "${value}"
    done 
}

section2() {
    values=($(seq 1 100))
    for value in "${values[@]}"; do 
        Rscript RESET_ensemble.R data sim2_"${value}".csv result_data/reset_all_sim2_"${value}".csv "${value}"
    done 
}

# Check command line arguments and execute the corresponding section

if [[ $# -eq 0 ]]; then
    echo "No section specified. Usage: RESET_ensemble.sh [section]"
    exit 1
fi

section=$1

case $section in
    section1)
        section1
        ;;
    section2)
        section2
        ;;
    *)
        echo "Invalid section. Available sections: section1, section2"
        exit 1
        ;;
esac
