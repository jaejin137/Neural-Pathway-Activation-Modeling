#!/bin/bash -l
#PBS -l walltime=10:00:00,pmem=2580mb,nodes=8:ppn=10
#PBS -m abe
#PBS -M slops001@umn.edu
#module load python2
#source activate myenv

#cd /panfs/roc/groups/2/johnsonm/slops001/Abbott_pulse_delays/Simulations/stim26/

#cd $HOME/Barb/Modeling_Data_rev1/Data_Processing/Barb_R_GPeSTN/Simulations

#module load intel # temporarily commented out by JL
python import_and_compile_ffem_multi.py
#mpirun -np 80 python simulate_parallel_MSI_ffem_multi.py
mpirun -np 32 python simulate_parallel_multi_ffem_multi.py --filename Neuron_Data # modified by JL
