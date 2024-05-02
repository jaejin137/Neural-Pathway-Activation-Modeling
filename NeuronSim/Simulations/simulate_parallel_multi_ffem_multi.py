#!/usr/bin/env python
"""Using the task-pull paradigm for high-throughput computingwith the 
module mpi4py. Task pull is an efficient way to perform a large number of
independent tasks when there are more tasks than processors, especially
when the run times vary for each task. 
This code is a modified version of code created by Craig Finch (cfinch@ieee.org).
Originally inspired by http://math.acadiau.ca/ACMMaC/Rmpi/index.html
"""

import os                    # to determine the number of neurons in your neuron_data folder
import numpy as np           # import np for arrays
from time import time        # timer for logging execution time
from NRTL_multi import Axon        # import Axon from NRTL meta class
from mpi4py import MPI       # use mpi4py for parallelization
from scipy.io import savemat # savemat for saving .mat files
from scipy.io import loadmat # loadmat for loading .mat files
from math import *

#-------------------------------------------------------------------------------
# Simulation settings:
# varibles here are global and therefore should not be expected to change during
# the simulation
#-------------------------------------------------------------------------------
#freq = 1000        # stimulation frequency (Hz)
delay = 20        # delay before and after stimulus train (msec): input waveform is zero during the delay
AP_delay = 20      # amount of time after the stimulation pulse during which APs will be looked for (msec)
nTest_pulses = 1 # number of stimulation pulses in the stimulus train

up_limit = 10   # maximum scaling factor and starting point for threshold search algorithm
epsilon = 0.001    # allowed error when searching for threshold
low_limit = 1e-3  # minimum scaling factor for threshold search algorithm

neuron_data_varible_name = "axon"      # variable name in each of your .mat files
#neuron_data_foldername = "filled later" # folder containing .mat files with neuron geom and voltage data
neuron_data_foldername = 'Neuron_Data' # folder containing .mat files with neuron geom and voltage data
geom_name=['GPeSTN', 'STNMCPath', 'ET035_R_Dent_2_VimE', 'ET035_R_Dent_2_VimI','ET035_R_Dent_2_Ret','ET035_R_SSC_2_VPL','ET035_R_pMC_2_VLA','ET035_R_ML' ,'ET035_R_MC_2_VimI','ET035_R_MC_2_VimE','ET035_R_LF','ET035_R_TF','ET035_R_IC','ET035_R_ZI']
#geom_name=['ET035_L_MC_2_VimE','ET035_L_MC_2_VimI','ET035_L_ML','ET035_L_pMC_2_VLA','ET035_L_SSC_2_VPL','ET035_L_Dent_2_Ret' ,'ET035_L_Dent_2_VimI','ET035_L_Dent_2_VimE','ET035_L_LF','ET035_L_TF','ET035_L_IC','ET035_L_ZI']

FFEM = False
Abbott_waveform = True # Originally, was "False". Chaged to "True" by JL. 
if FFEM: # Fourier FEM 
    dt = 0.001 # ms integration time-step for neuron and dt for waveform (msec) NOTE: your waveform will be downsampled as needed
    # IMPORTANT: Be sure your dt is no greater than 0.01 (i.e. Fs = 100 kHz) for 
    # voltage-regulated stimulation and no greater than 0.001 (i.e. Fs = 1 MHz) for
    # current-regulated stimulation (per Butson 2005). This is the sampling rate 
    # of the waveform and the integration time-step of the neuron simulation. 
    # If your simulations are using the Fourier FEM method use the code below to build
    # your waveform and set your dt. 
    Fs = 1024e3          # Sampling frequency [Hz] 
    #
    stim_amp = np.zeros((2,6)) # waveform stimulation amplitude (should always be one)
    stim_amp[0,:] = np.tile(-0.6667,6)
    stim_amp[1,:] = np.tile(-0.3333,6)

    if Abbott_waveform:
        freq=130
        # If your simulations are quasistatic and you are using a waveform stored 
        # as a .mat file use this code
        # filename of .mat with stimulation (should be located in NeuronModels/NeuronCode/waveforms
        #stim_pulse_filename = "./Abbott_recorded_waveform_short.mat" # waveform filename
        #stim_pulse_filename = "/home/jaejin/Barb/Modeling_Data_rev1/Data_Processing/Barb_R_GPeSTN/Simulations/Abbott_recorded_waveform.mat" # waveform filename
        stim_pulse_filename = "./Abbott_recorded_waveform.mat" # waveform filename
        stim_pulse_varible_name = "waveform" # name of the variable you'd like to load in filename.mat
        dt = 0.001                            # integration time-step for neuron and dt for waveform (msec)
        # IMPORTANT: Be sure your dt is no greater than 0.01 (i.e. Fs = 100 kHz). This 
        # is the sampling rate of the waveform and the integration time-step of the 
        # neuron simulation. 
        phase=0
        # load the stimulation pulse
        if isinstance(stim_pulse_filename,str):     
            temp = loadmat(stim_pulse_filename,squeeze_me=True,appendmat=True,struct_as_record=False) # load .mat file
            stim_pulse = np.tile(temp[stim_pulse_varible_name],(stim_amp.shape[0],1)) # load into an array the structure field: structure_field
            stim_pulse_base = stim_pulse[0,:] # load into an array the structure field: structure_field
    else: 
        # IMPORTANT NOTE: Your waveform will eventually get downsampled to 1000e3. Be sure this doesn't distort you waveform. You should be fine if starting with 1024e3.
        pw = round(60e6/Fs) # pulse width [us]
        # Build square pulse 
        # build waveform
        a = np.tile(0.0,pw)
        b = np.tile(-1,pw)
        c = np.tile(1,pw)
        d = np.tile(0.0,1024-(3*pw))
        stim_pulse_base = np.concatenate([a,b,c,d])
        #stim_pulse=np.concatenate([[stim_pulse],[stim_pulse],[stim_pulse]]) 
    
    delay_offset=np.zeros((2,6))
    delay_offset[0,:] = [0,500, 1000, 1500, 2500, 3750]
    delay_offset[1,:] = np.tile(0,6)

  
    freq=130
    #stim_pulse=np.zeros((phase.shape[0],stim_amp.shape[1],np.max(num_peaks)/(np.max(sine_freq)/1000)/dt+np.max(delay_offset)/dt+1)) # length of vector 1024 -> 1024 samples/1msec 
    #sine=np.zeros([phase.shape[0],5,stim_amp.shape[1]])      #peaks, frequency(kHz), phase offset (rad), delay(ms), stim_amp  
    stim_pulse=np.zeros((stim_amp.shape[0],int(len(stim_pulse_base)), stim_amp.shape[1]))
    for n in range(stim_amp.shape[1]): # trial number/amplitude combination -> 0:19
        for j in range(stim_amp.shape[0]): # electrode number -> 0:7
            stim_pulse[j,:,n]=stim_pulse_base*stim_amp[j,n]
       

        #sine[n,:,:]=np.concatenate([num_peaks,sine_freq,np.matrix(phase[n,:]),delay_offset[0,:],stim_amp],axis=0) ##peaks, freq (kHz), phase offset, delay offset, stim amplitude 
        # build waveform - contact 0 
        #scaled_stim_pulse_partial=np.zeros((stim_amp.shape[1],np.max(sine[n,0,:])/(np.max(sine[n,1,:])/1000)/dt+np.max(sine[n,3,:])/dt+1))
        #for i in range(stim_amp.shape[1]): 
        #   x = np.linspace(0,sine[n,0,i]*sine[n,1,i]/1000,sine[n,0,i]*(sine[n,1,i]/1000)/dt+1)   
        #   pulse = np.sin(2*pi*sine[n,1,i]/1000*x+sine[n,2,i])
        #   scaled_stim_pulse_partial[i,0:len(pulse)+sine[n,3,i]/dt]=np.concatenate([np.tile(0.0,sine[n,3,i]),pulse])*sine[n,4,i]
            
        #stim_pulse[n,:,:] = scaled_stim_pulse_partial
    
    stim_pulse_filename = False

else: # Quasistatic
    freq=130
    # If your simulations are quasistatic and you are using a waveform stored 
    # as a .mat file use this code
    
    # filename of .mat with stimulation (should be located in NeuronModels/NeuronCode/waveforms
    #stim_pulse_filename = "dbsCurrentControlledWaveform_136Hz_90usec_pseudoBiphasic_Lempka2010Recording" # waveform filename (commented out and replaced by JL)
    #stim_pulse_filename = "./Abbott_recorded_waveform_short.mat" # waveform filename
    #stim_pulse_filename = "/home/jaejin/Barb/Modeling_Data_rev1/Data_Processing/Barb_R_GPeSTN/Simulations/Abbott_recorded_waveform.mat" # waveform filename
    stim_pulse_filename = "./Abbott_recorded_waveform.mat" # waveform filename
    stim_pulse_varible_name = "waveform" # name of the variable you'd like to load in filename.mat
    dt = 0.01                            # integration time-step for neuron and dt for waveform (msec)
    # IMPORTANT: Be sure your dt is no greater than 0.01 (i.e. Fs = 100 kHz). This 
    # is the sampling rate of the waveform and the integration time-step of the 
    # neuron simulation. 
    
    
    #load stim_amp.mat file here#
    stim_amp_file = loadmat('stim_amp_all.mat',squeeze_me=True,appendmat=True,struct_as_record=False)
    stim_amp = stim_amp_file['stim_amp']
    
    
    
    phase=0
    # If the simulation is not Fourier FEM then load the stimulation pulse
    if type(stim_pulse_filename) is str:
        temp = loadmat(stim_pulse_filename,squeeze_me=True,appendmat=True,struct_as_record=False) # load .mat file
        stim_pulse = np.tile(temp[stim_pulse_varible_name],(stim_amp.shape[0],1)) # load into an array the structure field: structure_field   


#-------------------------------------------------------------------------------
# Enumerated for tags
#------------------------------------------------------------------------------- 
def enum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)
    
# Define MPI message tags
tags = enum('READY', 'DONE', 'EXIT', 'START')

# Initializations and preliminaries
comm = MPI.COMM_WORLD   # get MPI communicator object
size = comm.size        # total number of processes
rank = comm.rank        # rank of this process
status = MPI.Status()   # get MPI status object


#-------------------------------------------------------------------------------
# Simulation
#-------------------------------------------------------------------------------   
def simulate(neuron_filename,neuron_ind,stimulus,sine,stim_amp,delay_offset):#[neuron_list[iNeuron],iNeuron,n,stim_pulse,sine,stim_amp,delay_offset]
    
    # Run stimulation
    print("Neuron #" + str(neuron_ind+1))
    result={} # dictionary in which to store results
    
    # Load neuron data
    neuron_data_foldername = 'Neuron_Data' # folder containing .mat files with neuron geom and voltage data
    temp = loadmat('./' + neuron_data_foldername + '/' + neuron_filename,squeeze_me=True,appendmat=True,struct_as_record=False) # load .mat files
    neuron_data = temp[neuron_data_varible_name] # load into an array the structure field: structure_field
    
    # Build and initialize objects
    neuron = Axon() # create class instance
    neuron.construct_geom(neuron_data) # construct a multi-compartment cell or axon
    #neuron.initialize_simulation(dt,delay,freq,nTest_pulses,AP_delay,stimulus,neuron_data.voltage,sine,stim_amp,delay_offset)
    neuron.initialize_simulation(dt,delay,freq,nTest_pulses,AP_delay,stimulus,neuron_data.voltage,sine,stim_amp)
    
    print("Process keeps running even after error with neuron.initialization!!!")

    # Find stimulus threshold
    neuron.threshold(low_limit,up_limit,epsilon)
    result['thresh']=neuron.thresh
    result['nExtraAPs']=neuron.nExtraAPs
    
    # Find node where AP originated
    neuron.initial_AP_site() 
    result['time1']=neuron.time1
    result['time2']=neuron.time2
    result['AP_site1']=neuron.AP_site1
    result['AP_site2']=neuron.AP_site2
    
    # Calculate firing rate during stimulation
    neuron.calculate_firing_rate(neuron.delay,neuron.delay+neuron.dur) 
    result['avg_freq']=neuron.avg_freq
    result['avg_inst_freq']=neuron.avg_inst_freq
    result['std_inst_freq']=neuron.std_inst_freq
    
    # Stim wave and membrane voltage
    result['membrane_v']=np.array(neuron.membrane_v)
    result['membrane_v_time']=np.array(neuron.membrane_v_time)
    temp=np.array(neuron.v_ext)
    result['stim']=temp[neuron.recNode,:] ### Warning: neuron.recNode as an index to save the membrane voltage (neuron.v_ext) from the recording node. recNode is a node number, not a compartment number from allseclist. The indices from neuron.v_ext represent the sections in order as they were populated using allseclist. Therefore, what you actually need for the index is the index that represents the location of recNode in allseclist.
    result['stimulus_train']=np.array(neuron.stimulus_train) 
    result['scaled_stim_pulse_partial']=np.matrix(neuron.scaled_stim_pulse_partial) 
    result['scaled_stim_pulse']=np.array(neuron.scaled_stim_pulse) 
    result['stim_pulse']=np.array(neuron.stim_pulse) 
    result['scaled_electrodes']=np.array(neuron.electrodes) 
    #result['F3']=np.array(neuron.F3)
    #result['F4']=np.array(neuron.F4)
    #result['F5_3']=np.array(neuron.F5_3)
    
    # Finish simulation
    del neuron # deconstruct class instance
    return result # return the result dictionary

#-------------------------------------------------------------------------------
# Push-pull parallelization
#------------------------------------------------------------------------------- 
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--filename',action='store', dest="filename",type='str',help='what is the mat file that you want to run')
(opt,args) = parser.parse_args()
assert opt.filename, "filename not defined!"

neuron_data_foldername = opt.filename #"ET030_L_MC_2_VimI_200" # folder containing .mat files with neuron geom and voltage data
print(neuron_data_foldername) 

if rank == 0: # if rank is 0 make master
    # Master process executes code below

    
    # Get start time
    start = time()


    # make a list of .mat files containing geometry and voltage data for simulation
    neuron_list = sorted(os.listdir('./' + neuron_data_foldername)) # list of files containing geom and voltage info
    
    loop_end=stim_amp.shape[1]
    loop_start=0 #starting value of n
    n=loop_start    
    

    # Preallocate memory for list of results
    nNeuron = len(neuron_list) #Number of neurons
    all_results = [None]*(nNeuron*(loop_end-loop_start)) # list of results


    # If needed, create a folder for saving the results
    for j in range(loop_start,loop_end):   
        if not os.path.exists('results_'+ neuron_data_foldername + '_PSO' + str(j)):
            os.makedirs('results_'+ neuron_data_foldername + '_PSO' + str(j))
    #-------------------------------------------------------------------------------
    # Send and receive jobs to and from workers
    iNeuron = 0
    num_workers = size - 1
    closed_workers = 0
    print("Master starting with %d workers" % num_workers)
    
    while closed_workers < num_workers:
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status) #output, iNeuron, n 
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            # Worker is ready, so send it a task
            if iNeuron < nNeuron:
                if FFEM: 
                    comm.send([neuron_list[iNeuron],iNeuron,n,stim_pulse[:,:,n],0,np.array([1, 1]),delay_offset[:,n]], dest=source, tag=tags.START)#([neuron_list[iNeuron],iNeuron,n,stim_pulse[n,:,:]],sine[n,:,:],stim_amp,delay_offset],
                else: 
                    comm.send([neuron_list[iNeuron],iNeuron,n,stim_pulse,0,stim_amp[:,n],np.zeros((stim_amp.shape[0],1))], dest=source, tag=tags.START) #[neuron_list[iNeuron],iNeuron,n,stim_pulse,sine,stim_amp,delay_offset]
                print("Sending cell %d to worker %d" % (iNeuron, source))
                iNeuron += 1
                if (iNeuron==nNeuron and n<loop_end-1): 
                   n=n+1
                   iNeuron=0
                   print("iNeuron reset to %d" % iNeuron) 
            else:
                comm.send(None, dest=source, tag=tags.EXIT)
                print("Checkpoint 0") 
        elif tag == tags.DONE:
            all_results[data[1]+nNeuron*(data[2]-loop_start)] = data[0] #data[0]=output, data[1]=iNeuron, data[2]=n 
            savemat("results_"+ neuron_data_foldername + "_PSO" + str(data[2]) +"/" + str(neuron_list[data[1]]),all_results[data[1]+nNeuron*(data[2]-loop_start)],appendmat=True)
            print("Cell %d complete" % data[1])
        elif tag == tags.EXIT:
            print("Worker %d exited." % source)
            closed_workers += 1 
    print("While loop complete") 
    # Save results
    out_folder = os.getcwd().replace("Simulations", "Results")
    if not os.path.exists(out_folder):
        os.mkdir(out_folder) # create results folder
    for k in range(0,loop_end-loop_start):
        output=np.zeros((nNeuron,10)) # preallocate mem (nNeuron,10)
        for index in range(k*nNeuron,(k+1)*nNeuron):
            axon=neuron_list[index-k*nNeuron]
            iNeuron=int(axon[4:-4])
            output[iNeuron-1][0] = iNeuron
            output[iNeuron-1][1] = all_results[index]['thresh']
            output[iNeuron-1][2] = all_results[index]['nExtraAPs']
            output[iNeuron-1][3] = all_results[index]['avg_freq']
            output[iNeuron-1][4] = all_results[index]['avg_inst_freq']
            output[iNeuron-1][5] = all_results[index]['std_inst_freq']
            output[iNeuron-1][6] = all_results[index]['AP_site1']
            output[iNeuron-1][7] = all_results[index]['AP_site2']
            output[iNeuron-1][8] = all_results[index]['time1']
            output[iNeuron-1][9] = all_results[index]['time2']  
            # Save results    
            savemat(out_folder + "/Results_"+ neuron_data_foldername + "_PSO" + str(k+loop_start),{'output':output},appendmat=True)
        
        # Calculate and display elapsed time
    print("Elapsed time: " + str(time()-start) + " seconds")
    print("Master finishing")  

else: # if rank is not 0 make worker
    # Worker processes execute code below
    name = MPI.Get_processor_name()
    print("I am a worker with rank %d on %s." % (rank, name))
    while True:
        comm.send(None, dest=0, tag=tags.READY)
        task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)#([neuron_list[iNeuron],iNeuron,n,stim_pulse[n,:,:]],sine[n,:,:],stim_amp,delay_offset], dest=source, tag=tags.START)
        tag = status.Get_tag()
        if tag == tags.START:
            # Do the work here
            output=simulate(task[0],task[1],task[3],task[4],task[5],task[6]) #[neuron_list[iNeuron],iNeuron,n,stim_pulse,sine,stim_amp,delay_offset]
            comm.send([output,task[1], task[2]],dest=0, tag=tags.DONE)
        elif tag == tags.EXIT:
            break
    comm.send(None, dest=0, tag=tags.EXIT)

#delete neuron files to save space on MSI server 
#import shutil 
#shutil.rmtree(os.getcwd() + '/x86_64')
#import glob 
#for hgx in glob.glob(os.getcwd() +'/*.mod'):
#    os.remove(hgx)
