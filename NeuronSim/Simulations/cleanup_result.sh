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

#PWD="/home/jaejin/mnt/gdrive_jaejin/Work/PAM/Barb/Modeling_Data_rev1/Data_Processing/Barb_R_GPeSTN/Simulations"

if [ ! -d $PWD/$1 ]; then
    echo -n "Create a new directory, $1? (y/n)"
    read create_dir
    if [ "$create_dir" == "y" ]; then
        mkdir $PWD/$1;
    elif [ "$create_dir" == "n" ]; then
        echo "Directory, $1, already exists. Overwrite? (y/n)"
        read overwrite_dir
        if [ ! "$overwrite_dir" == "y" ]; then
            echo "Give me a new directory name. Exiting..."
            exit 1
        fi
    fi
fi

echo "Moving results to a new directory..."
mv $PWD/Neuron_Data $PWD/$1/
mv $PWD/*.png $PWD/$1/
#mv $PWD/for\ * $PWD/$1/
#mv $PWD/*.mod $PWD/$1/
#mv $PWD/*.inc $PWD/$1/
#mv $PWD/x86_64 *.png $PWD/$1/
#mv $PWD/__pycache__ $PWD/$1/
echo "Done"
echo
exit 0
