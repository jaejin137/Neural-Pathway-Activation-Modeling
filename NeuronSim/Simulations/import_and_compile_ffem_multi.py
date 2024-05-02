"""
Copy .mod files, NRTL metaclass, and waveform file (if applicable) 
to the current working directory and compile the NEURON .mod files.
"""  

import os       # for getting home and cwd directory
import sys      # for appending folders to path
import ctypes   # for message box
import platform # for getting the operating system
import py_compile           # for compiling NRTL
from shutil import copy     # for copying files
from subprocess import call # for running the command nrnivmodl
from distutils.dir_util import copy_tree # for copying folder contents

# User settings
stim_filename = True #False # Fourier FEM
#stim_filename = "dbsCurrentControlledWaveform_136Hz_90usec_pseudoBiphasic_Lempka2010Recording" # waveform filename
stim_filename = "Abbott_recorded_waveform" # Just for test (JL)


# Message box for windows users: instruct user to compile Neuron code manually
def win_compile_complete():
    cwd = os.getcwd() # current working directory (simulation directory)
    cmple_flag = []
    for iFile in os.listdir(cwd):
        if iFile.endswith('.mod'):
            file_name = os.path.splitext(iFile)[0]                
            cmple_flag.append(os.path.isfile(file_name+'.c'))

    if not cmple_flag:
        raise Exception("Error: No .mod files found")
        
    return cmple_flag

if __name__ == "__main__":
    
    if platform.system() == 'Windows': # if Windows OS
        
        # get repo location
        usr_name = os.environ.get( "USERNAME" )
        home = 'C:\\Users\\' + usr_name + '\\' # path to home
        #repo_loc = home + 'Documents\\GitHub\\NeuronModels\\NeuronCode' # GitHub repository location
        repo_loc = home + 'PAM\\NeuronModels\\NeuronCode' # GitHub repository location
        
        # get current working directory (simulation directory)
        cwd = os.getcwd() 
        
        # copy mod files
        mod_loc = repo_loc + '\\modFiles' # location of mod files
        copy_tree(mod_loc,cwd) # copy all mod files to current working directory
        
        # copy waveform to current directory if one is given
        if type(stim_filename) is str:
            wav_loc = repo_loc + '\\waveforms\\' + stim_filename + '.mat' # location of waveform files
            copy(wav_loc,cwd) # copy waveform to current working directory
        
        # import NRTL meta class
        NRTL_loc = repo_loc + '\\NRTLtest.py' # location of NRTL.py
        copy(NRTL_loc,cwd) # copy NRTL.py to current working directory
        
        # check if .mod files are compiled
        compile_flag = win_compile_complete()
        
        # instruct user to compile mechanisims if needed            
        while not all(compile_flag): 
            response = ctypes.windll.user32.MessageBoxA(0, "Sorry to see you're working on a windows machine. Your uncompiled .mod files have been copied to your simulation directory:\n\n" + os.getcwd() + "\n\nOpen mknrndll and compile the mod files in this folder then click OK","Windows :(", 1)
            if response==2:
                raise SystemExit(0)
            else:
                # check if .mod files are compiled
                compile_flag = win_compile_complete()

    else: # if Mac OS or Linux
        
        # get repo location
        home = os.path.expanduser('~') # path to home
        #repo_loc = home + '/GitHub/NeuronModels/NeuronCode' # GitHub repository location
        repo_loc = home + '/PAM/NeuronModels/NeuronCode' # GitHub repository location

        # get current working directory (simulation directory)
        cwd = os.getcwd() 

        # copy mod files
        mod_loc = repo_loc + '/modFiles' # location of mod files
        copy_tree(mod_loc,cwd) # copy all mod files to current working directory

        # copy waveform to current directory if one is given
        if type(stim_filename) is str:
            wav_loc = repo_loc + '/waveforms/' + stim_filename + '.mat' # location of waveform files
            copy(wav_loc,cwd) # copy waveform to current working directory

        # import NRTL meta class
        NRTL_loc = repo_loc + '/NRTL_multi.py' # location of NRTL.py
        #copy(NRTL_loc,cwd) # copy NRTL.py to current working directory
        py_compile.compile('NRTL_multi.py') # compile NRTL

        # compile mechanisims
        call(['nrnivmodl'])
        print('Finished compiling')
