#!/bin/bash

if [ $# -lt 1 ]; then
    echo
    echo "Target model name is missing. Exiting..."
    echo
    exit 1
fi

DIR_NEURONSIM=/home/jaejin/Data_Processing/NeuronSim

echo
echo "Chaging directory ..."
cd $DIR_NEURONSIM/$1/Analysis

echo
echo "Start activation analysis ..."
python analyze_population_jl.py

echo
echo "Back to NeuronSim directory"
cd $DIR_NEURONSIM
echo

