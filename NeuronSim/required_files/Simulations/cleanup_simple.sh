#!/bin/bash

# Check if all arguments are given
if [ $# -lt 1 ]; then
    echo
    echo "Usage: $0 [NEW DIR]"
    echo "  [NEW DIR]: a new directory to save all simulation results into."
    echo
    exit 1
fi

mkdir $1
mv *.png $1/
mv results_Neuron_Data_PSO* $1/
