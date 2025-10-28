#!/bin/bash

section1() {
    start=$1
    end=$2
    for value in $(seq "$start" "$end"); do
        Rscript AdaPT-GMM.R data sim1_"${value}".csv result_data/adapt_gmm_sim1_"${value}".csv "${value}"
    done
}

section2() {
    start=$1
    end=$2
    for value in $(seq "$start" "$end"); do
        Rscript AdaPT-GMM.R data sim2_"${value}".csv result_data/adapt_gmm_sim2_"${value}".csv "${value}"
    done
}

# Usage check
if [[ $# -lt 3 ]]; then
    echo "Usage: sh AdaPT-GMM.sh [section] [start] [end]"
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
