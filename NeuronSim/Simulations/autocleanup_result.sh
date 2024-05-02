#!/bin/bash
#
# Usage: cleanup_results [NEW DIR]
# [NEW DIR]: a new directory to save all simulation results into.

# Check if all arguments are given
if [ $# -lt 1 ]; then
    echo
    echo "Usage: $0 [NEW DIR]"
    echo "  [NEW DIR]: a new directory to save all simulation results into."
    echo
    exit 1
fi

if [ ! -d $1 ]; then
    mkdir $1;
fi

echo "Moving results to a new directory..."
mv *.png $1/
mv results_Neuron_Data_PSO* $1/
echo "Done"
echo
