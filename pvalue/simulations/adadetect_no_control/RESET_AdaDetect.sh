#!/bin/bash

section1() {
    values=($(seq 1 100))
    for value in "${values[@]}"; do 
        Rscript AdaDetect_RESET.R data1 sim1_"${value}".csv result_data1/reset_adadetect_sim1_"${value}".csv "${value}"
    done 
}

# Check command line arguments and execute the corresponding section

if [[ $# -eq 0 ]]; then
    echo "No section specified. Usage: RESET_AdaDetect.sh [section]"
    exit 1
fi

section=$1

case $section in
    section1)
        section1
        ;;
    *)
        echo "Invalid section. Available sections: section1"
        exit 1
        ;;
esac
