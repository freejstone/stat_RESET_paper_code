#!/bin/bash

section1() {
    values=($(seq 1 50))
    for value in "${values[@]}"; do 
        Rscript AdaPT.R data1 sim1_"${value}".csv gam result_data1/adapt_sim1_"${value}".csv "${value}" 's(x1, x2)' 's(x1, x2)'
    done 
}

section2() {
    values=($(seq 51 100))
    for value in "${values[@]}"; do 
        Rscript AdaPT.R data1 sim1_"${value}".csv gam result_data1/adapt_sim1_"${value}".csv "${value}" 's(x1, x2)' 's(x1, x2)'
    done 
}

section3() {
    values=($(seq 1 50))
    for value in "${values[@]}"; do 
        Rscript AdaPT.R data1 sim2_"${value}".csv gam result_data1/adapt_sim2_"${value}".csv "${value}" 's(x1, x2)' 's(x1, x2)'
    done 
}

section4() {
    values=($(seq 51 100))
    for value in "${values[@]}"; do 
        Rscript AdaPT.R data1 sim2_"${value}".csv gam result_data1/adapt_sim2_"${value}".csv "${value}" 's(x1, x2)' 's(x1, x2)'
    done 
}

section5() {
    values=($(seq 1 50))
    for value in "${values[@]}"; do 
        Rscript AdaPT.R data1 sim3_"${value}".csv gam result_data1/adapt_sim3_"${value}".csv "${value}" 's(x1, x2)' 's(x1, x2)'
    done 
}

section6() {
    values=($(seq 51 100))
    for value in "${values[@]}"; do 
        Rscript AdaPT.R data1 sim3_"${value}".csv gam result_data1/adapt_sim3_"${value}".csv "${value}" 's(x1, x2)' 's(x1, x2)'
    done 
}

section7() {
    values=($(seq 1 50))
    for value in "${values[@]}"; do 
        Rscript AdaPT.R data2 sim_"${value}".csv glmnet result_data2/adapt_sim_"${value}".csv "${value}"
    done 
}

section8() {
    values=($(seq 51 100))
    for value in "${values[@]}"; do 
        Rscript AdaPT.R data2 sim_"${value}".csv glmnet result_data2/adapt_sim_"${value}".csv "${value}"
    done 
}

# Check command line arguments and execute the corresponding section

if [[ $# -eq 0 ]]; then
    echo "No section specified. Usage: AdaPT.sh [section]"
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
    section5)
        section5
        ;;
    section6)
        section6
        ;;
    section7)
        section7
        ;;
    section8)
        section8
        ;;
    *)
        echo "Invalid section. Available sections: section1, section2, section3, section4, section5, section6, section7, section8"
        exit 1
        ;;
esac
