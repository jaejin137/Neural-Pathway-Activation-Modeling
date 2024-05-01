#!/bin/bash

if [ $# -lt 1 ]; then
    echo
    echo "Target model name is missing. Exiting..."
    echo
    exit 1
fi

DIR_AXONS=/home/jaejin/Data_Processing/Axons/Interp_$1/Barb_R_GPeSTN/16axons
DIR_NEURONSIM=/home/jaejin/Data_Processing/NeuronSim

echo
echo "Starting simulation ..."

echo
echo "Creating simulation directory ..."
mkdir $DIR_NEURONSIM/$1

echo
echo "Copying requried files ..."
cp -r $DIR_NEURONSIM/required_files/* $DIR_NEURONSIM/$1/

echo
echo "Copying axon files ..."
cp -r $DIR_AXONS $DIR_NEURONSIM/$1/Simulations/Neuron_Data

echo
echo "Chaging directory ..."
cd $DIR_NEURONSIM/$1/Simulations

echo
echo "Running Neuron simulation ..."
sh RunNeuron.sh

echo
echo "Simulation completed."

echo
echo "Cleanig up ..."
sh autocleanup_result.sh results_n16
echo
echo "Clean-up completed."
echo
echo "Back to NeuronSim directory"
cd $DIR_NEURONSIM
echo

