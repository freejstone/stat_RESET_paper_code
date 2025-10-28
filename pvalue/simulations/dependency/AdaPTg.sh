#!/bin/bash

section1() {
    start=$1
    end=$2
    for value in $(seq "$start" "$end"); do
        Rscript AdaPTg.R data sim1_"${value}".csv gam result_data/adaptg_sim1_"${value}".csv "${value}" 's(x)' 's(x)'
    done
}

section2() {
    start=$1
    end=$2
    for value in $(seq "$start" "$end"); do
        Rscript AdaPTg.R data sim2_"${value}".csv gam result_data/adaptg_sim2_"${value}".csv "${value}" 's(x)' 's(x)'
    done
}

# Usage check
if [[ $# -lt 3 ]]; then
    echo "Usage: sh AdaPTg.sh [section] [start] [end]"
    exit 1
fi

section=$1
start=$2
end=$3

case $section in
    section1)
        section1 "$start" "$end"
        ;;
    section2)
        section2 "$start" "$end"
        ;;
    *)
        echo "Invalid section. Available sections: section1, section2"
        exit 1
        ;;
esac
