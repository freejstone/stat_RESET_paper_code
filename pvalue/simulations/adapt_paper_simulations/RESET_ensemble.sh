#!/bin/bash

section1() {
    values=($(seq 1 1))
    for value in "${values[@]}"; do 
        Rscript RESET_ensemble.R data1 sim1_"${value}".csv result_data1/reset_sim1_CHECK"${value}".csv "${value}"
    done 
}

section2() {
    values=($(seq 1 1))
    for value in "${values[@]}"; do 
        Rscript RESET_ensemble.R data1 sim2_"${value}".csv result_data1/reset_sim2_CHECK"${value}".csv "${value}"
    done 
}

section3() {
    values=($(seq 1 1))
    for value in "${values[@]}"; do 
        Rscript RESET_ensemble.R data1 sim3_"${value}".csv result_data1/reset_sim3_CHECK"${value}".csv "${value}"
    done 
}

section4() {
    values=($(seq 1 1))
    for value in "${values[@]}"; do 
        Rscript RESET_ensemble.R data2 sim_"${value}".csv result_data2/reset_sim_CHECK"${value}".csv "${value}"
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
    section3)
        section3
        ;;
    section4)
        section4
        ;;
    *)
        echo "Invalid section. Available sections: section1, section2, section3, section4"
        exit 1
        ;;
esac
