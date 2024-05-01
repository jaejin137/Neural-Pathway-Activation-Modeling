from math import pi,sqrt    # import some math functions
from neuron import h        # import neuron class
import numpy as np          # for working with arrays
import scipy.signal as sig # for downsampling when running Fourier FEM
import os                    # os for getting current working directory
import platform             # for getting the operating system
from abc import ABCMeta
from cmath import exp 
from scipy.io import savemat # savemat for saving .mat files

# Load mechanisms manually if running on a windows OS (otherwise mechs are loaded automatically)
if platform.system() == 'Windows':
    h.nrn_load_dll(os.getcwd() + '\\nrnmech.dll')    # load mechanisims

# if Linux OS assume user is on the supercomputer and prevent need for X11
if platform.system() == 'Linux':  
    import matplotlib
    matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt # plotting functions
else:
    import pylab as plt # import plotting functions normally


class NRTL():
    """
    The NRTL metaclass contains methods that are inherited by each of the child
    classes (e.g. Axon,GPi_Cell).
       
    Methods:     
        initialize_simulation - passes simulation relevant attributes from a StimulusTrain object
        insert_vext - calculated and plays the extracellular voltages into the multicompartment neuron
        pulse_train - create a train of pulses from user inputs: pulse and frequency
        integrate - solve the multicompartment neuron equivalent circuit model at each time point using NERUON 
        response_to_stimulus - determine if a neuron was activated
        trial - insert_vext, integrate, and response_to_stimulus in that order
        threshold - determine the neuron threshold for activation by running multiple trials
        initial_AP_site - identify the node where from an action potiential originated
        calculate_firing_rate - calculated the firing rate of a firing axon
        
    NEURON Globals:
        h.dt:                 (double) integration time-step for neuron (msec)
    
    Attributes:
        FEM_in_mV:            (numpy array of doubles) voltage at each compartment. Will be 2D (nCompartments x nFrequencies) if using Fourier FEM
        delay:                (int or double) delay before and after stimulus train (msec): input waveform is zero during the delay
        freq:                 (int or double)stimulation frequency (Hz): the number of samples in your input waveform should be int(1e5/freq)-1
        nTest_pulses:         (int)number of stimulation pulses in the stimulus train
        dur:                  (double) duration for which stim pulses are delivered (msec)
        tstop:                (double) time for the entire stimulus train (msec)
        AP_delay:             (int or double) amount of time after the stimulation pulse during which APs will be looked for (ms)
        stim_pulse_begin:     (list of doubles) stim pulse times
        nIdeal_APs:           (int) ideal number of APs (one for each stim pulse) 
        nReq_APs:             (int) required APs to be considered Activated (80%)
        stimulus_train:       (list of doubles) stimulus pulse train (volts)
        stimulus_train_time:  (list of doubles) time points for stimulus pulse train (msec)
        thresh:               (double) minimum stim amplitude needed for activation 
    """

    __metaclass__ = ABCMeta
    
    
    def initialize_simulation(self,dt,delay,freq,nTest_pulses,AP_delay,stim_pulse,FEM_in_V,sine=0,electrodes=[1, 0, 1, 0]):

        """
        Pass values to the cell or axon object, calculate windows where we'll be
        looking to find action potientials, and import FEM results i.e. 
        compartment voltages.
        
        >>> cell_obj.initialize_simulation(dt,delay,freq,nTest_pulses,AP_delay,stim_pulse,FEM_in_V,sine,electrodes=0)
        
        Inputs:
            dt:           (double) integration time-step for neuron and dt for waveform (msec)
            delay:        (int or double) delay before and after stimulus train (msec): input waveform is zero during the delay
            freq:         (int or double) stimulation frequency (Hz): the number of samples in your input waveform should be int(1e5/freq)-1
            nTest_pulses: (int) number of stimulation pulses in the stimulus train
            AP_delay:     (int or double) amount of time after the stimulation pulse during which APs will be looked for (msec)
            make_plot:    (bool) True to generate a plot of the stimulus train. Default is False
        
        Attributes:
            dt:           (double) integration time-step for neuron (msec)
            delay:        (int or double) delay before and after stimulus train (msec): input waveform is zero during the delay
            freq:         (int or double)stimulation frequency (Hz): the number of samples in your input waveform should be int(1e5/freq)-1
            nTest_pulses: (int)number of stimulation pulses in the stimulus train
            dur:          (double) duration for which stim pulses are delivered (msec)
            tstop:        (double) time for the entire stimulus train (msec)
            AP_delay:     (int or double) amount of time after the stimulation pulse during which APs will be looked for (ms)
            nIdeal_APs:   (int) ideal number of APs (one for each stim pulse) 
            nReq_APs:     (int) required APs to be considered Activated (80%)
            stim_pulse_begin:    (list of doubles) stim pulse times
            stimulus_train:      (list of doubles) stimulus pulse train (volts)
            stimulus_train_time: (list of doubles) time points for stimulus pulse train (msec)
            sine: variables for building a sinusoid waveform: nx5 matrix (n = number of electrodes), [#peaks, frequency(Hz), phase offset (rad), delay(ms), amplitude] 
	    electrodes: indicates which electrodes are turned on 1 = on, 0 = off
        """
            
        # Setup extracellular stimulation constants
        h.dt = dt                           # time step for integration (also time step for waveform)
        self.delay = delay                  # delay before and following stimulus train
        self.freq = freq                    # frequency of stim pulse delivery
        self.nTest_pulses = nTest_pulses    # number of stim pulses delivered
        self.dur = nTest_pulses*1000.0/freq # duration for which stim pulses are delivered (ms)
        self.tstop = delay+self.dur+delay   # time for the entire stimulus train (ms) 
        self.AP_delay = AP_delay            # amount of time after the stimulation pulse during which APs will be looked for (ms)
        self.stim_pulse = stim_pulse        # import stimulation pulse 
        self.FEM_in_V = FEM_in_V            # import FEM results
        self.sine=sine 			    # variables for building a sinusoid waveform: nx5 matrix (n = number of electrodes), [#peaks, frequency(Hz), phase offset (rad), delay(ms), amplitude]		
        self.electrodes=electrodes          # matrix indicating which electrodes are turned on 1 = on (and polarity ex. -1 or +1), 0 = off
        
              
        #-----------------------------------------------------------------------
        # Calculate the starting time of each stimulus pulse, and store in the
        # Vector required_APs. APs should follow after a certain time whose 
        # maximum is AP_delay.
        #-----------------------------------------------------------------------
        # initialize list of stim pulse times
        self.stim_pulse_begin = []
        
        # calculate the number of required APs and make sure it matches the number of test pulses  
        if self.freq == 1:
            self.nIdeal_APs = 1 # if the frequency is one, look for only one AP
        else:
            if self.dur > self.tstop:
                self.nIdeal_APs = int((self.tstop-self.delay)*self.freq/1000)  # if the duration is larger than tstop then subtract out the delay
            else:
                self.nIdeal_APs = int(self.dur*self.freq/1000)        
    
        for iIdeal_APs in range(self.nIdeal_APs):
            self.stim_pulse_begin.append(self.delay+iIdeal_APs*1000/self.freq)
        
        # number of APs required for Activation is 80% of the number of stim pulses
        self.nReq_APs = int(self.nIdeal_APs*1.0)
        
        print("      Stimulation frequency = " , str(self.freq) , " Hz" )
        print("      # of stimulation pulses = " , str(self.nTest_pulses)) 
        print("      delay during which APs will be recorded = " , str(self.AP_delay))
        print("      ideal # of APs = " , str(self.nIdeal_APs))
        print("      # of APs required for Activation = " , str(self.nReq_APs))
        
        # error checking 
        if self.nIdeal_APs != self.nTest_pulses:
            print( "Number of APs is not " , str(self.nTest_pulses) + "!!, but " , str(self.nIdeal_APs))
            
        
        # Plot stimulus pulse train
        len(stim_pulse) 
        plt.figure()  
        if len(stim_pulse.shape)==2: #multiple-electrodes
            stimulus_train=np.zeros([stim_pulse.shape[0],stim_pulse.shape[1]])            
            for i in range(stim_pulse.shape[0]): 
                #stimulus_train=self.pulse_train(stim_pulse[i,:])*self.electrodes[0,i] # original
                stimulus_train=self.pulse_train(stim_pulse[i,:])*self.electrodes[i] # by JL for debug
                stimulus_train_time = list(ii*h.dt for ii in range(0,stimulus_train.shape[0])) # simulation time points
                plt.plot(stimulus_train_time,stimulus_train)
                plt.xlabel('Time (ms)')
                plt.ylabel('Amplitude')
                plt.title('Stimulation train')
                plt.grid(True)
                name="stimulation_waveform_%i.png" %i
                plt.savefig(name, format='png')  
                plt.clf()
                
        else:  #single electrode 
            stimulus_train=self.pulse_train(stim_pulse)
            stimulus_train_time = list(ii*h.dt for ii in range(0,stimulus_train.shape)) # simulation time points
            plt.plot(stimulus_train_time,stimulus_train)
            plt.xlabel('Time (ms)')
            plt.ylabel('Amplitude')
            plt.title('Stimulation train')
            plt.grid(True)
            plt.savefig("stimulation_waveform.png")  
            plt.close("all")
            
        
    #---------------------------------------------------------------------------
    # Calculate extracellular voltage and play into cell
    #---------------------------------------------------------------------------   
    def insert_vext(self,stim_amplitude): 
        """
        Calculate extracellular voltage for each compartment at each time point
        and play them into the e_extracellular mechanism of their respective 
        compartments.
        
        >>> cell_obj.initialize_simulation(stim_amplitude)
        
        Inputs:
            stim_amplitude: (double) stimulation scaling factor for simulation
        """  
                                                              
        # create some list of extracellular potentials on each segment and time vector
        print ("                 Trial scaling factor = " , str(stim_amplitude))

        v_ext = []
        iCompartment = 0
        # play v_ext into e_extracellular reference


        if (len(self.FEM_in_V.shape) == 1 and isinstance(self.sine,np.ndarray)): #fourier, sinusoid, single-electrode, single frequency 
            print( "Fourier, single electrode, sine")                    
        elif (len(self.FEM_in_V.shape) == 2 and self.FEM_in_V.shape[1] > 0 and isinstance(self.sine,np.ndarray)): # Fourier, sine, multi electrode (#electrodes x 1 row x compartments)                   
            print( "Fourier, multi electrode,sine")
        elif (len(self.FEM_in_V.shape) == 1 and self.FEM_in_V.shape[1] == 1 and isinstance(self.sine,int)): # non- Fourier, single electrode (1 electrode x 1 row x compartments) 
            print( "non-Fourier, single electrode")
        elif (len(self.FEM_in_V.shape) == 2 and self.FEM_in_V.shape[1] > 0 and isinstance(self.sine,int)): # non-Fourier, multi electrode (#electrodes x 1 row x compartments) 
            print("non-Fourier, multi electrode" )

        for sec in self.allseclist:
            # if the array of voltages is 1D then scale the waveform by FEM_in_V 
            # that was calculated at each compartment
            if (len(self.FEM_in_V.shape) == 1 and isinstance(self.sine,np.ndarray)): #fourier, sinusoid, single-electrode, single frequency 
                scaled_stim_pulse_partial=np.zeros((1,self.sine[0]/(self.sine[1]/1000)/h.dt+self.sine[3]/h.dt+1))
                time = np.linspace(0,self.sine[0]*self.sine[1]/1000,self.sine[0]*self.sine[1]/1000/h.dt+1)   
                pulse = abs(self.FEM_in_V[iCompartment])*np.sin(2*pi*self.sine[1]/1000*time+self.sine[2]+np.angle(self.FEM_in_V[iCompartment]))
                scaled_stim_pulse=np.concatenate([np.tile(0.0,self.sine[3]),pulse])*stim_amplitude*1000
	        #print "Fourier, single electrode, sine"                    
            elif (len(self.FEM_in_V.shape) == 2 and self.FEM_in_V.shape[1] > 0 and isinstance(self.sine,np.ndarray)): # Fourier, sine, multi electrode (#electrodes x 1 row x compartments)                   
                scaled_stim_pulse_partial=np.zeros((self.FEM_in_V.shape[1],np.max(self.sine[0,:])/(np.max(self.sine[1,:])/1000)/h.dt+np.max(self.sine[3,:])/h.dt+1))
                for i in range(self.FEM_in_V.shape[1]): 
                    time = np.linspace(0,self.sine[0,i]*self.sine[1,i]/1000,self.sine[0,i]*self.sine[1,i]/1000/h.dt+1)   
                    pulse = abs(self.FEM_in_V[iCompartment,i])*np.sin(2*pi*self.sine[1,i]/1000*time+self.sine[2,i]+np.angle(self.FEM_in_V[iCompartment,i]))
                    scaled_stim_pulse_partial[i,self.sine[3,i]/h.dt:len(pulse)+self.sine[3,i]/h.dt]=np.concatenate([np.tile(0.0,self.sine[3,i]),pulse])*self.electrodes[i]
                scaled_stim_pulse = sum(scaled_stim_pulse_partial,0)*stim_amplitude*1000
	        #print "Fourier, multi electrode,sine"
                
            elif (len(self.FEM_in_V.shape) == 1 and self.FEM_in_V.shape[1] == 1 and isinstance(self.sine,int)): # non- Fourier, single electrode (1 electrode x 1 row x compartments) 
               # scale pulse by stim_amplitude and convert FEM result from V to mV
               scaled_stim_pulse = self.stim_pulse* self.FEM_in_V[iCompartment]*stim_amplitude*1000
	       #print "non-Fourier, single electrode" 
               
            elif (len(self.FEM_in_V.shape) == 2 and self.FEM_in_V.shape[1] > 0 and isinstance(self.sine,int)): # non-Fourier, multi electrode (#electrodes x 1 row x compartments) 
               scaled_stim_pulse_partial=np.zeros((self.FEM_in_V.shape[1],self.stim_pulse.shape[1]))
               for i in range(self.FEM_in_V.shape[1]): 
                   scaled_stim_pulse_partial[i,:] = self.stim_pulse[i,:]* self.FEM_in_V[iCompartment,i]*self.electrodes[i]
               scaled_stim_pulse = sum(scaled_stim_pulse_partial,0)*stim_amplitude*1000 
	       #print "non-Fourier, multi electrode" 
	       #print "stim" + str(self.electrodes) 
                
            # if the array of voltages is 2D then perform Fourier FEM
            elif (len(self.FEM_in_V.shape) == 2 and isinstance(self.sine,int)): # fourier, single-electrode ( #frequencies x #compartments x 1 electrode )         
                                  
                scaled_raw_pulse = self.stim_pulse * stim_amplitude  # scale pulse by stim_amplitude

                # Fourier FEM
                n = 1024;                            # 1024 point DFT
                L = self.stim_pulse.shape[1]             # Length of signal
                Y= np.fft.fft(scaled_raw_pulse,n);   # perform fft
            
                # reflect FEM results
                F1 = np.concatenate([self.FEM_in_V[:,iCompartment], np.conj(self.FEM_in_V[np.arange((n/2)-1,0,-1),iCompartment])]) 
                F2 = F1*Y                     # scale and phase shift fft by FEM result
                F3 = np.fft.ifft(F2,n)        # perform inverse fft
                F4 = np.real(F3)              # remove small imaginary aspects that remain due to computational precision
                F5 = sig.resample(F4,int(1/h.dt))  # resample to match h.dt
                
                scaled_stim_pulse = F5*1000 # convert result from V to mV and perform superposition               
                                                           
            elif (len(self.FEM_in_V.shape) == 3 and self.FEM_in_V.shape[0] > 1 and isinstance(self.sine,int)): # fourier, multi-electrode ( #frequencies x #compartments x 1 electrode )                                                           
                # scale pulse by stim_amplitude
		#print "Fourier, multi electrode" 
                F5=np.zeros((self.stim_pulse.shape[0],1000))
                for i in range(self.FEM_in_V.shape[2]):
                    scaled_raw_pulse = self.stim_pulse[i,:] * stim_amplitude

                    # Fourier FEM
                    n = 1024;                            # 1024 point DFT
                    L = self.stim_pulse.shape[1]         # Length of signal
                    Y= np.fft.fft(scaled_raw_pulse,n);   # perform fft
                
                    # reflect FEM results
                    F1 = np.concatenate([self.FEM_in_V[:,iCompartment,i], np.conj(self.FEM_in_V[np.arange((n/2)-1,0,-1),iCompartment,i])]) 
                    F2 = F1*Y                     # scale and phase shift fft by FEM result
                    F3 = np.fft.ifft(F2,n)        # perform inverse fft
                    F4 = np.real(F3)              # remove small imaginary aspects that remain due to computational precision
                    F5[i,:] = sig.resample(F4,int(1/h.dt))  # resample to match h.dt
                    
                scaled_stim_pulse =sum(F5,0)*1000  # convert result from V to mV and perform superposition 
                scaled_stim_pulse_partial= F5
            
                    
                                                                                                                                                                                                                                                                                      
            else:
                raise Exception("Voltage value array has more than two dimensions")   
            
            # construct a stimulus pulse train from a scaled_stim_pulse
            self.stimulus_train = self.pulse_train(scaled_stim_pulse)
            self.scaled_stim_pulse=scaled_stim_pulse
            self.scaled_stim_pulse_partial=scaled_stim_pulse_partial 
            #self.F3=F3 
            #self.F4=F4
            #self.F5=F5
                
            v_ext.append(h.Vector(self.stimulus_train))
            iCompartment += 1 # for indexing into the voltage
            

            for seg in sec:
                v_ext[-1].play(seg._ref_e_extracellular, h.dt)
        self.v_ext = v_ext 
          
    #---------------------------------------------------------------------------
    # Construct a pulse train from a single pulse at self.freq Hz
    #---------------------------------------------------------------------------           
    def pulse_train(self,stimulation_pulse):    
        """
        Create a train of stimulation pulses from the user provided pulse and 
        stimulation frequency
        
        >>> cell_obj.pulse_train(stimulation_pulse)
        
        Inputs:
            stimulation_pulse: (numpy array of doubles) a single stimulation pulse, normally has zeros on either side
        """  
        # initialize stimulus pulse train at zero
        nStimulus_train_time = int(self.tstop/h.dt) + 1 # number of simulation time points
        stimulus_train = [0] * nStimulus_train_time
    
        # insert waveform to stimulus train
        my_time = self.delay                         # delay in msec
        my_indx = int(self.delay/h.dt) + 1           # index at delay time point
        nCycle_size = int(1000/(self.freq*h.dt)) + 1 # number of time points in a single pulse given the frequency and
         
        # starting from "delay" add one  
        while (my_time < self.delay+self.dur and my_indx < nStimulus_train_time): # during stimulation time
            for iCycle_size in range(nCycle_size):
                if iCycle_size < len(stimulation_pulse): # if iCycle_size is equal to or larger then the recorded waveform does not contain voltage values for the points (this will happen if you are stimulating at a lower frequency than the waveform was recorded at)
                    stimulus_train[my_indx] = stimulation_pulse[iCycle_size]  # set stimulus for simulation

                my_indx += 1 
                my_time += h.dt
                
        return np.array(stimulus_train)
           
           
    #---------------------------------------------------------------------------
    # Solve equations
    #---------------------------------------------------------------------------
    def integrate(self):
        """
        Solve the multicompartment neuron equivalent circuit model at each time
        point using NEURON.
        
        >>> cell_obj.initialize_simulation()
        """       
        h.finitialize(self.v_init)
        h.fcurrent()
        while h.t < self.tstop:
            h.fadvance()



    #---------------------------------------------------------------------------
    # Checks if cell responds to all stimulus pulses and records any extra APs 
    # to a vector extra_APs
    #---------------------------------------------------------------------------
    def response_to_stimulus(self):
        """
        Determine if a neuron was activated.
        
        >>> cell_obj.initialize_simulation()
        """  
        nAPs = 0                 # number of required APs that happened
        nExtra = 0                 # number of required APs that happened
        self.extra_AP_times = [] # number of extra APs that happened
        flag = 0 
        flag2 = 0

        for ii in range(int(self.apc.n)):
            # AP_delay added if there are any APs that are responding to stimulus 
            # pulse at the very end of pulse duration
            if (self.apc_times.x[ii] > self.delay) and (self.apc_times.x[ii] <= self.delay+self.dur+self.AP_delay):
                flag = 0  
                for jj in range(int(flag2),self.nIdeal_APs):
                    if (self.apc_times.x[ii] > self.stim_pulse_begin[jj]) and (self.apc_times.x[ii] <= self.stim_pulse_begin[jj]+self.AP_delay):
                        nAPs += 1
                        flag = 1
                        flag2 = jj + 1
                        break
                        
                if flag == 0:
                    self.extra_AP_times.append(self.apc_times.x[ii])
                    nExtra += 1
    
        # length of extra_AP_times list is the number of extra APs
        self.nExtraAPs = nExtra #len(self.extra_AP_times) 

        print( "      # of stimulation pulses = " , str(self.nTest_pulses)) 
        print( "      # of required APs = " , str(self.nReq_APs) )
        print( "      # of APs recorded within " , str(self.AP_delay) , " msec post-stimulus window = " + str(nAPs))
        print( "      # of extra APs recorded within " , str(self.AP_delay) , " msec post-stimulus window = " + str(self.nExtraAPs))
    
        if (nAPs >= self.nReq_APs): # there is AP following 80% of stimulus pulses
            print( "      AP following at least 80% of stimulus pulses? Yes")
            return True
            
        print( "      AP following at least 80% of stimulus pulses? No")
        return False
            

    def trial(self,stim_amplitude):
        """
        Run a trial: Calculate and play extracellular voltage into each 
        solve equivalent circuit model, and check the cell response.
        
        >>> cell_obj.initialize_simulation(stim_amplitude)
        
        Inputs:
            stim_amplitude: (double) stimulation scaling factor for simulation
        """      
        self.insert_vext(stim_amplitude)
        self.integrate()
        return self.response_to_stimulus() # returns 1 if the number of APs needed was recorded       

    #-------------------------------------------------------------------------------
    # Return voltage threshold for stimulation
    #-------------------------------------------------------------------------------
    def threshold(self,low_limit,up_limit,epsilon):
        """
        Find minimum threshold for Activation using a binary search algorithm 
        and running multiple simulations. 
        
        >>> cell_obj.threshold(low_limit,up_limit,epsilon)
        
        Inputs:
            low_limit: (double) minimum scaling factor
            up_limit:  (double) maximum scaling factor
            epsilon:   (double) allowed error when searching for threshold
        """       
        lbound = low_limit
        ubound = up_limit
        strength = ubound
        
        while abs(ubound-lbound) > abs(epsilon*ubound):
            excited = self.trial(strength)
            
            if excited:
                ubound = strength
                strength = (ubound+lbound)/2
            else:
                if strength == up_limit:
                    self.thresh = -1 # -1 if the axon did not fire at the upper limit
                    return
                else:
                    lbound = strength
                    strength = (3*ubound+lbound)/4
        
        if not excited:
            self.trial(ubound) # run one more time at the max amplitude that evokes an AP to get other measurements
        
        # return strength  # threshold in Volts
        # strength is between ubound and lbound; ubound is returned because that is the 
        # closest value to strength for which we know causes an AP (strength might not be enough)  
        #return ubound
        self.thresh = ubound
        return
        
        
    #---------------------------------------------------------------------------
    # Find AP initiation site
    #---------------------------------------------------------------------------
    def initial_AP_site(self):
        """
        Find AP initiation site.
        
        >>> cell_obj.initial_AP_site()
        
        """   

        # Initialize variables
        self.AP_site1 = -1   	#node where AP appears first
        self.AP_site2 = -1   	#node where AP appears second
        self.time1=-1
        self.time2=-1
        firstAP_node_number = []
        firstAP_times = []
               
        # For each AP counter get the index and time of occurrence for the first AP that was recorded after the delay
        for iAPcounters in range(self.nAPcounters):     # for each APcounter (there is one on each node)
            AP_times = self.AP_timerec[iAPcounters]     # retrieve times that each recorded AP occurred from AP_counters
            nAPs_total = len(AP_times)                  # total number of APs that occurred

            # For each recorded AP on a given APcounter (iAPcounters) retrieve the time of the first AP then break out of the loop
            for iAPs_total in range(nAPs_total):                   # for each recorded AP
                if AP_times[iAPs_total] > self.delay:            # if an AP was recorded after the delay
                    firstAP_node_number.append(iAPcounters)        # store node number
                    firstAP_times.append(AP_times[iAPs_total])   # retrieve the time of the first AP on that node
                    break                                         # and break out of the loop
        
        # if the firstAP_nodes list is still empty, notify user that there were no recorded APs
        if not firstAP_times:
            print( "NO AP during stimulation....")
            
        # otherwise get the node number and time for the 1st and second APs
        else:
            # sort firstAP_times and firstAP_node_number from lowest time value to highest
            sorted_firstAP_node_number = [x for (y,x) in sorted(zip(firstAP_times,firstAP_node_number))]
            sorted_firstAP_times = sorted(firstAP_times)
            
            # find index of AP counter that recorded shortest time (i.e. where AP first appeared)
            self.time1 = sorted_firstAP_times[0]
            self.AP_site1 = sorted_firstAP_node_number[0]
            
            if len(sorted_firstAP_times)>1: # if the cell fires more than 1 AP
                # find second site for AP initialization
                self.time2 = sorted_firstAP_times[1]
                self.AP_site2 = sorted_firstAP_node_number[1]
    

        
    #---------------------------------------------------------------------------
    # Calculate the average firing rate, period, and variance of each for the 
    # axon or cell during a given time block 
    #---------------------------------------------------------------------------
    def calculate_firing_rate(self, time_beg, time_end):
        """
        Calculate the average firing rate, period, and variance of each for the 
        axon or cell during a given time block.
        
        >>> cell_obj.calculate_firing_rate(time_beg, time_end)
        
        Inputs:
            time_beg: (double) beginning time (ms) 
            time_end: (double) end time (ms)
        """   
        ISI_sum = 0
        ISI_squared_sum = 0
        inst_freq_sum = 0
        inst_squared_freq_sum = 0
        nISI = 0 # number of inter-spike intervals
        
        nRecorded_APs = int(self.apc.n) # number of recorded APs
        for iRecorded_APs in range(nRecorded_APs-1): # for each recorded AP
            # if the recorded AP is between time_beg and time_end
            if (self.apc_times.x[iRecorded_APs] > time_beg) and (self.apc_times.x[iRecorded_APs+1] <= time_end):
                ISI = self.apc_times.x[iRecorded_APs+1]-self.apc_times.x[iRecorded_APs] # calculate ISI
                ISI_sum += ISI
                ISI_squared_sum += ISI**2
                inst_freq_sum += (1000/ISI) # convert from ms to seconds, then to Hz
                inst_squared_freq_sum += (1000/ISI)**2
                nISI+=1
            
        self.avg_period = 0
        self.avg_inst_freq = 0
        if nISI != 0:
            self.avg_period = ISI_sum/nISI
            self.avg_inst_freq = inst_freq_sum/nISI
            
        nSpikes = 0 # number of spikes
        for iRecorded_APs in range(nRecorded_APs):
            if (self.apc_times.x[iRecorded_APs] > time_beg) and (self.apc_times.x[iRecorded_APs] <= time_end):
                nSpikes+=1
            
        self.avg_freq = 1000*nSpikes/(time_end-time_beg) # convert time from ms to sec
        self.std_period = 0
        self.std_inst_freq = 0
        if nISI-1 > 0:
            var_period = (ISI_squared_sum - ISI_sum**2/nISI)/(nISI-1)
            var_freq = (inst_squared_freq_sum - inst_freq_sum**2/nISI)/(nISI-1)
            if var_period > 0:
                self.std_period = sqrt(var_period)
            if var_freq > 0:
                self.std_inst_freq = sqrt(var_freq)
        

class STN(NRTL):
    """
    This STN class inherits the NRTL metaclass methods.
    
    Methods:     
    construct_geom - build cell geometry (also calls create_seclist and insert_APcounters
    create_seclist - creates a list of all the sections
    insert_APcounters - insert action potential counters and membrane voltage recorders
    
    Parent (NRTL) Methods:     
    simulation initialization - passes simulation relevant attributes from a StimulusTrain object
    insert_vext - calculated and plays the extracellular voltages into the multicompartment neuron
    integrate - solve the multicompartment neuron equivalent circuit model at each time point using NEURON 
    response_to_stimulus - determine if a neuron was Activated
    trial - insert_vext, integrate, and response_to_stimulus in that order
    threshold - determine the neuron threshold for Activation by running multipule trials
    
    NEURON Globals:
    h.dt:           (double) integration time-step for neuron (msec)
    h.celsius:      (double) temperature in Celsius
    
    Attributes:
    v_init:          (double) cell initial voltage
    nStinPerStretch: (int) number of STIN in a row
    nSoma:           (int) number of compartments of type "soma"
    nDend:           (int) number of compartments of type "dend"
    nInitseg:        (int) number of compartments of type "initseg"
    nNodes:          (int) number of comparmtents of type "node"
    nMysa:           (int) number of comparmtents of type "mysa"
    nFlut:           (int) number of comparmtents of type "flut"
    nStin:           (int) number of comparmtents of type "stin"
    soma:            (NEURON Section object) one for each "soma" compartment
    dend:            (NEURON Section object) one for each "dend" compartment
    initseg:         (NEURON Section object) one for each "initseg" compartments
    node:            (NEURON Section object) one for each "node" compartment
    MYSA:            (NEURON Section object) one for each "mysa" compartment
    FLUT:            (NEURON Section object) one for each "flut" compartment
    STIN:            (NEURON Section object) one for each "stin" compartment
    allseclist       (NERUON list) of sections
    apc:             (NEURON APCount object) for counter APs
    apc_times:       (NEURON Vector object) for storing times of APs
    nAPcounters:     (int) number of AP counters, equal to nNodes          
    AP_counters:     (list of NEURON APCount objects) one on each node
    AP_timerec:      (NEURON Vector object) for storing times of APs
    membrane_v:      (NEURON Vector object) for recording the membrane voltage on the last node
    """
    
    def construct_geom(self,STNstrip,fiberD=2):        
        """
        Build a multicompartment cell, generate a section list, and insert AP 
        counters 
        
        >>> axon_obj.construct_geom(fiberD,STNstrip)
        
        Inputs:
            fiberD:    (double) fiber diameter (um): must be 1 or 2
            STNstrip: (structure) imported from a .mat and contains all the compartment information
        """
        # Set axon parameters and build axons
        h.celsius = 36  # temperature in Celsius
        self.v_init = -60 

        # get axon compartment data
        soma = STNstrip.soma
        dend = STNstrip.dend
        initseg = STNstrip.initseg
        node = STNstrip.node
        mysa = STNstrip.mysa
        flut = STNstrip.flut
        stin = STNstrip.stin
        #self.nStinPerStretch = STNstrip.nStinPerStretch # this was in laura's PPN code, but we don't need it for STN
        
        # calculate number of each compartment type
        self.nSoma = int(soma.shape[0]) # number of soma
        self.nDend = int(dend.shape[0]) # number of dendrites
        #self.nInitseg = int(initseg.shape[0]) # number of initial segments
        try: # check if only one initial segment
            initseg.shape[1]
        except IndexError:
            x=None
        
        if x is None:
            self.nInitseg = 1
        else:
            self.nInitseg = int(initseg.shape[0])
            
        self.nNodes = int(node.shape[0]) # number of nodes
        self.nMysa = int(mysa.shape[0]) # number of mysa
        self.nFlut = int(flut.shape[0]) # number of flut
        self.nStin = int(stin.shape[0]) # number of internodes
     
        # Axon electrical parameters
        rhoa=0.7e6 # Ohm-um
        mycm=0.1   # uF/cm2/lamella membrane
        mygm=0.001 # S/cm2/lamella membrane
            
        #-----------------------------------------------------------------------
        # Set morphological parameters based on fiber diameter
        #-----------------------------------------------------------------------
        space_p1 = 0.002  
        space_p2 = 0.004
        space_i = 0.004
        if fiberD==1: axonD=0.8; nodeD=0.7; paraD1=0.7; paraD2=0.8; deltax=200; nl=20
        elif fiberD==2: axonD=1.6; nodeD=1.4; paraD1=1.4; paraD2=1.6; deltax=200; nl=30
        else: raise Exception("Unknown fiber diameter")
        
        Rpn0=(rhoa*.01)/(pi*((((nodeD/2)+space_p1)**2)-((nodeD/2)**2)))
        Rpn1=(rhoa*.01)/(pi*((((paraD1/2)+space_p1)**2)-((paraD1/2)**2)))
        Rpn2=(rhoa*.01)/(pi*((((paraD2/2)+space_p2)**2)-((paraD2/2)**2)))
        Rpx=(rhoa*.01)/(pi*((((axonD/2)+space_i)**2)-((axonD/2)**2)))
        
        #-----------------------------------------------------------------------
        # Create "sections" for each compartment and assign mechanisms
        #-----------------------------------------------------------------------
        
        # Intracellular ion concentrations are typical mammalian values (from Johnson & Wu, 1999 via NEURON tutorial)
        cai0_ca_ion = 1e-4
        ki0_k_ion = 140
        nai0_na_ion = 10

        #Extracellular ion concentrations taken from Bevan,Wilson 1999 paper (slice bathing solution)	
        cao0_ca_ion = 2
        ko0_k_ion = 2.5
        nao0_na_ion = 126
        
        # parameters (mho/cm2) 
        my_gcaT_CaT = 0
        my_gcaL_HVA = 9.5e-4 
        my_gcaN_HVA = 1.15e-3
        my_gk_sKCa = 1.8*6.84e-5 
       	my_gk_KDR = 3.84e-3  
        my_gk_Kv31 = 1.2*1.34e-2
       	my_gna_Na = 0.75*1.48e-2   
       	my_gna_NaL =0.75*1.11e-5  
       	my_gk_Ih = 1.01e-3
       	my_gpas_STh = 7.84112e-05 
       	
       	# create self.nSoma soma compartments
       	self.soma = [h.Section() for ii in range(self.nSoma)]
 
       	
       	for ii in range(0,int(self.nSoma)):
            self.soma[ii].nseg = 1
            self.soma[ii].Ra = 150.2
            self.soma[ii].cm = 1
            # Passive Membrane properties
            #self.soma[ii].insert('myions') # Dummy mechanism to set up ion concentrations for Na and K ########################## TODO: laura you took this out. 
            # Autonomous fast firing membrane properties
            self.soma[ii].insert('cacum') # Model for CA2+ Accumulation 
            h.cai0_cacum = cai0_ca_ion
            self.soma[ii].insert('CaT') # Model for CaT
            self.soma[ii].insert('HVA') # Model for HVA Ca2+ Current
            self.soma[ii].insert('sKCa')
            self.soma[ii].insert('KDR') # Model for K+ Delayed Rectifier Current (fast-deactivating)
            self.soma[ii].insert('Kv31') 
            self.soma[ii].insert('Na') # Model for NA+ Current
            self.soma[ii].insert('NaL') # Model for NA+ Leak Current  
            self.soma[ii].insert('Ih') # 
            self.soma[ii].insert('STh')
            self.soma[ii].insert('extracellular')   
            
            for seg in self.soma[ii]:
                seg.gcaT_CaT = my_gcaT_CaT
                seg.HVA.gcaL = my_gcaL_HVA
                seg.HVA.gcaN = my_gcaN_HVA
                seg.gk_sKCa = my_gk_sKCa
                seg.gk_KDR = my_gk_KDR
                seg.gk_Kv31 = my_gk_Kv31
                seg.gna_Na = my_gna_Na
                seg.gna_NaL = my_gna_NaL                 
                seg.gk_Ih = my_gk_Ih
                seg.gpas_STh = my_gpas_STh
                seg.xraxial[0] = 1e09
                seg.xg[0] = 1e09
                seg.xc[0] = 0


                     
       	# create self.nDend dendrite compartments
       	self.dend = [h.Section() for ii in range(self.nDend)]
       	for ii in range(0,int(self.nDend)):
            self.dend[ii].nseg = 1
            self.dend[ii].Ra = 150.2
            self.dend[ii].cm = 1
            # Passive Membrane properties
            # self.dend[ii].insert('myions') # Dummy mechanism to set up ion concentrations for Na and K
            # Autonomous fast firing membrane properties
            self.dend[ii].insert('cacum') # Model for CA2+ Accumulation
            h.cai0_cacum = cai0_ca_ion
            self.dend[ii].insert('CaT') # Model for CaT
            self.dend[ii].insert('HVA') # Model for HVA Ca2+ Current
            self.dend[ii].insert('sKCa')
            self.dend[ii].insert('KDR') # Model for K+ Delayed Rectifier Current (fast-deactivating)
            self.dend[ii].insert('Kv31') 
            self.dend[ii].insert('Na') # Model for NA+ Current
            self.dend[ii].insert('NaL') # Model for NA+ Leak Current  
            self.dend[ii].insert('Ih') # 
            self.dend[ii].insert('STh') 
            self.dend[ii].insert('extracellular')            

            for seg in self.dend[ii]:
                seg.gna_Na = 0.65*1e-7
                seg.gna_NaL = 0.65*8.1e-6                 
                seg.gpas_STh = my_gpas_STh
                seg.xraxial[0] = 1e09 ########## seg.xraxial, seg.xg and seg.xc were not set for dend for some reason!
                seg.xg[0] = 1e09
                seg.xc[0] = 0
                
                
       	# create self.nInitseg initial segment compartment
       	self.initseg = [h.Section() for ii in range(self.nInitseg)]
       	
       	for ii in range(0,int(self.nInitseg)):
            self.initseg[ii].nseg = 1
            self.initseg[ii].Ra = 150.2
            self.initseg[ii].cm = 1
            # Passive Membrane properties
            # self.initseg[ii].insert('myions') # Dummy mechanism to set up ion concentrations for Na and K
            # Autonomous fast firing membrane properties
            self.initseg[ii].insert('cacum') # Model for CA2+ Accumulation
            h.cai0_cacum = cai0_ca_ion
            self.initseg[ii].insert('CaT') # Model for CaT
            self.initseg[ii].insert('HVA') # Model for HVA Ca2+ Current
            self.initseg[ii].insert('sKCa')
            self.initseg[ii].insert('KDR') # Model for K+ Delayed Rectifier Current (fast-deactivating)
            self.initseg[ii].insert('Kv31') 
            self.initseg[ii].insert('Na') # Model for NA+ Current
            self.initseg[ii].insert('NaL') # Model for NA+ Leak Current  
            self.initseg[ii].insert('Ih') # 
            self.initseg[ii].insert('STh')
            self.initseg[ii].insert('extracellular')            
            
            for seg in self.initseg[ii]:
                seg.gcaT_CaT = my_gcaT_CaT
                seg.HVA.gcaL = my_gcaL_HVA
                seg.HVA.gcaN = my_gcaN_HVA
                seg.gk_sKCa = my_gk_sKCa
                seg.gk_KDR = my_gk_KDR
                seg.gk_Kv31 = my_gk_Kv31
                seg.gna_Na = my_gna_Na
                seg.gna_NaL = my_gna_NaL                 
                seg.gk_Ih = my_gk_Ih
                seg.gpas_STh = my_gpas_STh
                seg.xraxial[0] = 1e09
                seg.xg[0] = 1e09
                seg.xc[0] = 0


        # create self.nNodes node compartments
        self.node = [h.Section() for ii in range(self.nNodes)] 
        
        # create all nodes that are not last node and insert mechanisms
        for ii in range(0,int(self.nNodes)-1):
            self.node[ii].nseg = 1
            self.node[ii].Ra = rhoa/10000
            self.node[ii].cm = 2
            self.node[ii].insert('axnode80') # a new mod file with correct STN values
            self.node[ii].insert('extracellular')
            #seg.gnabar_axnode75 = 2.0   
            #seg.gnapbar_axnode75 = 0.05
            #seg.gkbar_axnode75 = 0.07
            #seg.gl_axnode75 = 0.005
            #seg.ek_axnode75 = -85 
            #seg.ena_axnode75 = 55 
            #seg.el_axnode75 = -60
            #seg.vshift_axnode75 = 15 
            #seg.vtraub_axnode75 = -80 
            for seg in self.node[ii]:
                seg.xraxial[0] = Rpn0 # MOhms/cm
                seg.xg[0] = 1e10
                seg.xc[0] = 0
        
        # create last node (passive) and insert mechanisms  
        self.node[self.nNodes-1].nseg = 1 # number of segments
        self.node[self.nNodes-1].Ra = rhoa/10000 # axial resistivity (Ohm-cm)
        self.node[self.nNodes-1].cm = 2
        self.node[self.nNodes-1].insert('pas') # insert passive current
        self.node[self.nNodes-1].insert('extracellular') # insert extracellular mechanism
        
        for seg in self.node[self.nNodes-1]:
            seg.pas.g = 0.0001 # S/cm^2
            seg.pas.e = -65 
            seg.xraxial[0] = Rpn0 # MOhms/cm
            seg.xg[0] = 1e10
            seg.xc[0] = 0
        
        
        # create self.nMysa mysa compartments and insert mechanisms
        self.MYSA = [h.Section() for ii in range(0,self.nMysa)] 
        
        for ii in range(0,int(self.nMysa)):
            self.MYSA[ii].nseg = 1
            self.MYSA[ii].Ra = rhoa/10000
            self.MYSA[ii].cm = 2
            self.MYSA[ii].insert('pas')
            self.MYSA[ii].insert('extracellular') # insert extracellular mechanism
            for seg in self.MYSA[ii]:
                seg.pas.g = 0.0001 # S/cm^2
                seg.pas.e = -65 
                seg.xraxial[0] = Rpn1 # MOhms/cm
                seg.xg[0] = mygm/(nl*2)
                seg.xc[0] = mycm/(nl*2)	
        
        
        # create self.nFlut flut compartments and insert mechanisms
        self.FLUT = [h.Section() for ii in range(0,self.nFlut)]  
        
        for ii in range(0,int(self.nFlut)):
            self.FLUT[ii].nseg = 1
            self.FLUT[ii].Ra = rhoa/10000
            self.FLUT[ii].cm = 2
            self.FLUT[ii].insert('parak80')
            self.FLUT[ii].insert('pas') # insert passive current
            self.FLUT[ii].insert('extracellular') # insert extracellular mechanism
            
            for seg in self.FLUT[ii]:
                #seg.gkbar_parak75 = 0.02 
                #seg.ek_parak75 = -85
                #seg.vshift_parak75 = 15
                seg.pas.g = 0.0001 # S/cm^2
                seg.pas.e = self.v_init # mV
                seg.xraxial[0] = Rpn2 # MOhms/cm
                seg.xg[0] = mygm/(nl*2)
                seg.xc[0] = mycm/(nl*2)
        
        
        # create self.nStin stin compartments and insert mechanisms
        self.STIN = [h.Section() for ii in range(0,self.nStin)]  
        for ii in range(0,int(self.nStin)):
            self.STIN[ii].nseg = 1
            self.STIN[ii].Ra = rhoa/10000
            self.STIN[ii].cm = 2
            self.STIN[ii].insert('pas') # insert passive current
            self.STIN[ii].insert('extracellular') # insert extracellular mechanism
            
            for seg in self.STIN[ii]:
                seg.pas.g = 0.0001 # S/cm^2
                seg.pas.e = -65 
                seg.xraxial[0] = Rpx # MOhms/cm
                seg.xg[0] = mygm/(nl*2)
                seg.xc[0] = mycm/(nl*2)

        
        
        #---------------------------------------------------------------------------
        # Connect the sections in order to make a single multi-compartment cell with soma, dendrites, initial segment, and axons
        #   axon order: mysa, flut, stin, stin, stin, flut, mysa 
        #---------------------------------------------------------------------------
        # connect segments
        # hoc: connect child(cx),parent(px)
        # python: child.connect(parent,px,cx)
        
        # Somas
        for ii in range(0,2):
            self.soma[ii+1].connect(self.soma[ii],(1),(0))
            
        self.dend[0].connect(self.soma[0],(0),(0))
        
        #Dendrites
        for ii in range(0,18):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
        
        self.dend[19].connect(self.dend[1],(1),(0))
        
        for ii in range(19,32):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[33].connect(self.dend[19],(1),(0))    
        
        for ii in range(33,39):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[40].connect(self.dend[0],(1),(0))   
        
        for ii in range(40,49):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[50].connect(self.dend[40],(1),(0))   
        
        for ii in range(50,61):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[62].connect(self.dend[50],(1),(0))   
        
        for ii in range(62,69):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[70].connect(self.dend[68],(1),(0))  
        
        for ii in range(70,76):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[77].connect(self.soma[2],(1),(0))   
            
        for ii in range(77,86):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[87].connect(self.dend[77],(1),(0))   
        
        for ii in range(87,98):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[99].connect(self.dend[87],(1),(0))   
        
        for ii in range(99,112):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
                      
        self.dend[113].connect(self.dend[105],(1),(0))
         
        for ii in range(113,119):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[120].connect(self.soma[2],(1),(0))   
        
        for ii in range(120,131):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[132].connect(self.dend[124],(1),(0))   
         
        for ii in range(132,140):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[141].connect(self.dend[121],(1),(0))   
        
        for ii in range(141,155):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[156].connect(self.dend[120],(1),(0))   
        
        for ii in range(156,166):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[167].connect(self.soma[2],(1),(0))   
        
        for ii in range(167,176):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[177].connect(self.dend[173],(1),(0))   
        
        for ii in range(177,179):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[180].connect(self.dend[172],(1),(0))
        
        self.dend[181].connect(self.dend[169],(1),(0)) 
        
        for ii in range(181,189):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[190].connect(self.dend[168],(1),(0))   
        
        for ii in range(190,194):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[195].connect(self.dend[167],(1),(0))

        for ii in range(195,201):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[202].connect(self.dend[198],(1),(0))   

        for ii in range(202,206):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[207].connect(self.dend[197],(1),(0))   

        for ii in range(207,215):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[216].connect(self.dend[196],(1),(0))   

        for ii in range(216,225):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[226].connect(self.dend[220],(1),(0))   

        for ii in range(226,232):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[233].connect(self.dend[195],(1),(0)) 

        for ii in range(233,239):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[240].connect(self.soma[2],(1),(0))   

        for ii in range(240,249):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[250].connect(self.dend[244],(1),(0))   

        for ii in range(250,253):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[254].connect(self.dend[250],(1),(0))   

        for ii in range(254,260):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[261].connect(self.soma[0],(0),(0))   

        for ii in range(261,271):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[272].connect(self.soma[0],(0),(0))   

        for ii in range(272,276):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[277].connect(self.soma[0],(0),(0))  

        for ii in range(277,282):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[283].connect(self.dend[281],(1),(0))   
       
        self.dend[284].connect(self.dend[280],(1),(0))  

        self.dend[285].connect(self.dend[277],(1),(0))
        
        for ii in range(285,287):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))        

        self.initseg[0].connect(self.soma[2],(1),(0))
        self.node[0].connect(self.initseg[0],(1),(0)) 
        
        #Axons
        #STN axon
        self.MYSA[0].connect(self.node[0],(1),(0))
        self.FLUT[0].connect(self.MYSA[0],(1),(0))
        self.STIN[0].connect(self.FLUT[0],(1),(0))
        self.STIN[1].connect(self.STIN[0],(1),(0))
        self.STIN[2].connect(self.STIN[1],(1),(0))
        self.FLUT[1].connect(self.STIN[2],(1),(0))
        self.MYSA[1].connect(self.FLUT[1],(1),(0))
        self.node[1].connect(self.MYSA[1],(1),(0))
        
        for ii in range(1,self.nNodes-1): 
            self.MYSA[2*ii].connect(self.node[ii],(1),(0))
            self.FLUT[2*ii].connect(self.MYSA[2*ii],(1),(0))
            self.STIN[3*ii].connect(self.FLUT[2*ii],(1),(0))
            self.STIN[3*ii+1].connect(self.STIN[3*ii],(1),(0))
            self.STIN[3*ii+2].connect(self.STIN[3*ii+1],(1),(0))
            self.FLUT[2*ii+1].connect(self.STIN[3*ii+2],(1),(0))
            self.MYSA[2*ii+1].connect(self.FLUT[2*ii+1],(1),(0))
            self.node[ii+1].connect(self.MYSA[2*ii+1],(1),(0))

        
        #---------------------------------------------------------------------------
        # Create the actual geometry by placing each compartment at its correct 3-D 
        # coordinate 
        #---------------------------------------------------------------------------
        
        for ii in range(self.nSoma):
            h.pt3dadd(soma[ii,0],soma[ii,1],soma[ii,2],soma[ii,3],sec=self.soma[ii])
            h.pt3dadd(soma[ii,4],soma[ii,5],soma[ii,6],soma[ii,7],sec=self.soma[ii])
        for ii in range(self.nDend):
            h.pt3dadd(dend[ii,0],dend[ii,1],dend[ii,2],dend[ii,3],sec=self.dend[ii])
            h.pt3dadd(dend[ii,4],dend[ii,5],dend[ii,6],dend[ii,7],sec=self.dend[ii])
        if self.nInitseg ==1:
            h.pt3dadd(initseg[0],initseg[1],initseg[2],initseg[3],sec=self.initseg[0])
            h.pt3dadd(initseg[4],initseg[5],initseg[6],initseg[7],sec=self.initseg[0])
        else:
            for ii in range(self.nInitseg):
                h.pt3dadd(initseg[ii,0],initseg[ii,1],initseg[ii,2],initseg[ii,3],sec=self.initseg[ii])
                h.pt3dadd(initseg[ii,4],initseg[ii,5],initseg[ii,6],initseg[ii,7],sec=self.initseg[ii])
        for ii in range(self.nNodes):
            h.pt3dadd(node[ii,0],node[ii,1],node[ii,2],node[ii,3],sec=self.node[ii])
            h.pt3dadd(node[ii,4],node[ii,5],node[ii,6],node[ii,7],sec=self.node[ii])
        for ii in range(self.nMysa):
            h.pt3dadd(mysa[ii,0],mysa[ii,1],mysa[ii,2],mysa[ii,3],sec=self.MYSA[ii])
            h.pt3dadd(mysa[ii,4],mysa[ii,5],mysa[ii,6],mysa[ii,4],sec=self.MYSA[ii])
        for ii in range(self.nFlut):
            h.pt3dadd(flut[ii,0],flut[ii,1],flut[ii,2],mysa[ii,3],sec=self.FLUT[ii])
            h.pt3dadd(flut[ii,4],flut[ii,5],flut[ii,6],mysa[ii,7],sec=self.FLUT[ii])
        for ii in range(self.nStin):
            h.pt3dadd(stin[ii,0],stin[ii,1],stin[ii,2],stin[ii,3],sec=self.STIN[ii])
            h.pt3dadd(stin[ii,4],stin[ii,5],stin[ii,6],stin[ii,7],sec=self.STIN[ii])
        
        
        # Create a list of all the sections
        self.allseclist = self.create_seclist()
        
        # Create action potential (AP) counters and membrane voltage recorders
        self.insert_APcounters()
       


    #---------------------------------------------------------------------------
    # Create a list of all the sections
    #---------------------------------------------------------------------------
    def create_seclist(self):
        allseclist = h.SectionList()
        for ii in range(self.nSoma):
            allseclist.append(sec=self.soma[ii])
        for ii in range(self.nDend):
            allseclist.append(sec=self.dend[ii])         
        if self.nInitseg == 1:
            allseclist.append(sec=self.initseg[0])
        else:
            for ii in range(self.nInitseg):
                allseclist.append(sec=self.initseg[ii])
        for ii in range(self.nNodes):
            allseclist.append(sec=self.node[ii])
        for ii in range(self.nMysa):
            allseclist.append(sec=self.MYSA[ii])
        for ii in range(self.nFlut):
            allseclist.append(sec=self.FLUT[ii])
        for ii in range(self.nStin):
            allseclist.append(sec=self.STIN[ii]) 
        return allseclist
            
    #---------------------------------------------------------------------------
    # Create action potential (AP) counters and membrane voltage recorders
    #---------------------------------------------------------------------------
    def insert_APcounters(self):
        self.recNode = 10#######################DEBUG: we use node 11 (in python is 10) here because of debugging only. This value *should* be 27, which means we are recording from node 28, which is the furthest away from soma
        # Create a single action potential counter with a "thresh" setting in order 
        # to capture action potentials and determine the threshold for the axon
        self.apc = h.APCount(0.5,sec=self.node[self.recNode])
        self.apc_times = h.Vector()     # empty vector for action potential times
        self.apc.thresh = -20           # threshold that must be crossed to count as an action potential (mV)
        self.apc.record(self.apc_times) # set it to record

        # Set up AP counters in every node and branches (if axon has any)
        self.nAPcounters = self.nNodes
        
        # records AP time at every node
        self.AP_timerec = []
        for ii in range(self.nAPcounters): 
            self.AP_timerec.append(h.Vector())
    
        # Setup AP counter
        self.AP_counters = []
        for ii in range(self.nNodes):
            self.AP_counters.append(h.APCount(0.5,sec=self.node[ii]))    # put AP counter in every node
            self.AP_counters[ii].record(self.AP_timerec[ii])                  # records times of APs
        
        # Setup membrane voltage recorder at each time step
        self.membrane_v = h.Vector()
        self.membrane_v_time = h.Vector()
        self.membrane_v.record(self.node[self.recNode](0.5)._ref_v)
        self.membrane_v_time.record(h._ref_t)
        

class PPNType2(NRTL):
    """
    This PPNType2 class inherits the NRTL metaclass methods.
    
    Methods:     
    construct_geom - build cell geometry (also calls create_seclist and insert_APcounters
    create_seclist - creates a list of all the sections
    insert_APcounters - insert action potential counters and membrane voltage recorders
    
    Parent (NRTL) Methods:     
    simulation initialization - passes simulation relevant attributes from a StimulusTrain object
    insert_vext - calculated and plays the extracellular voltages into the multicompartment neuron
    integrate - solve the multicompartment neuron equivalent circuit model at each time point using NEURON 
    response_to_stimulus - determine if a neuron was Activated
    trial - insert_vext, integrate, and response_to_stimulus in that order
    threshold - determine the neuron threshold for Activation by running multiple trials
    
    NEURON Globals:
    h.dt:           (double) integration time-step for neuron (msec)
    h.celsius:      (double) temperature in Celsius
    
    Attributes:
    v_init:          (double) cell initial voltage
    nStinPerStretch: (int) number of STIN in a row
    nSoma:           (int) number of compartments of type "soma"
    nDend:           (int) number of compartments of type "dend"
    nInitseg:        (int) number of compartments of type "initseg"
    nNodes:          (int) number of comparmtents of type "node"
    nMysa:           (int) number of comparmtents of type "mysa"
    nFlut:           (int) number of comparmtents of type "flut"
    nStin:           (int) number of comparmtents of type "stin"
    soma:            (NEURON Section object) one for each "soma" compartment
    dend:            (NEURON Section object) one for each "dend" compartment
    initseg:         (NEURON Section object) one for each "initseg" compartments
    node:            (NEURON Section object) one for each "node" compartment
    MYSA:            (NEURON Section object) one for each "mysa" compartment
    FLUT:            (NEURON Section object) one for each "flut" compartment
    STIN:            (NEURON Section object) one for each "stin" compartment
    allseclist       (NERUON list) of sections
    apc:             (NEURON APCount object) for counter APs
    apc_times:       (NEURON Vector object) for storing times of APs
    nAPcounters:     (int) number of AP counters, equal to nNodes          
    AP_counters:     (list of NEURON APCount objects) one on each node
    AP_timerec:      (NEURON Vector object) for storing times of APs
    membrane_v:      (NEURON Vector object) for recording the membrane voltage on the last node
    """
    
    def construct_geom(self,cell_data,fiberD=2):        
        """
        Build a multicompartment cell, generate a section list, and insert AP 
        counteres 
        
        >>> axon_obj.construct_geom(fiberD,cell_data)
        
        Inputs:
            fiberD:    (double) fiber diameter (um): must be 1 or 2
            cell_data: (structure) imported from a .mat and contains all the compartment information
        """   
        # Set axon parameters and build axons
        h.celsius = 36  # temperature in Celsius
        self.v_init = -70  # axon resting potential

        # get axon compartment data
        soma = cell_data.soma
        dend = cell_data.dend
        initseg = cell_data.initseg
        node = cell_data.node
        mysa = cell_data.mysa
        flut = cell_data.flut
        stin = cell_data.stin
        self.nStinPerStretch = cell_data.nStinPerStretch
        
        # calculate number of each compartment type
        self.nSoma = int(soma.shape[0]) # number of soma
        self.nDend = int(dend.shape[0]) # number of dendrites
        #self.nInitseg = int(initseg.shape[0]) # number of initial segments
        
        try: # check if only one initial segment
            initseg.shape[1]
        except IndexError:
            x=None
        
        if x is None:
            self.nInitseg = 1
        else:
            self.nInitseg = int(initseg.shape[0])
            
        self.nNodes = int(node.shape[0]) # number of nodes
        self.nMysa = int(mysa.shape[0]) # number of mysa
        self.nFlut = int(flut.shape[0]) # number of flut
        self.nStin = int(stin.shape[0]) # number of internodes
     
        # Axon electrical parameters
        rhoa=0.7e6 # Ohm-um
        mycm=0.1   # uF/cm2/lamella membrane
        mygm=0.001 # S/cm2/lamella membrane
            
        #-----------------------------------------------------------------------
        # Set morphological parameters based on fiber diameter
        #-----------------------------------------------------------------------
        space_p1 = 0.002  
        space_p2 = 0.004
        space_i = 0.004
        if fiberD==1: axonD=0.8; nodeD=0.7; paraD1=0.7; paraD2=0.8; deltax=200; nl=20
        elif fiberD==2: axonD=1.6; nodeD=1.4; paraD1=1.4; paraD2=1.6; deltax=200; nl=30
        else: raise Exception("Unknown fiber diameter")
        
        Rpn0=(rhoa*.01)/(pi*((((nodeD/2)+space_p1)**2)-((nodeD/2)**2)))
        Rpn1=(rhoa*.01)/(pi*((((paraD1/2)+space_p1)**2)-((paraD1/2)**2)))
        Rpn2=(rhoa*.01)/(pi*((((paraD2/2)+space_p2)**2)-((paraD2/2)**2)))
        Rpx=(rhoa*.01)/(pi*((((axonD/2)+space_i)**2)-((axonD/2)**2)))
        
        #-----------------------------------------------------------------------
        # Create "sections" for each compartment and assign mechanisms
        #-----------------------------------------------------------------------
        
        # Intracellular ion concentrations are typical mammalian values (from Johnson & Wu, 1999 via NEURON tutorial)
        cai0_ca_ion = 1e-4
        ki0_k_ion = 140
        nai0_na_ion = 10

        #Extracellular ion concentrations taken from Nakanishi 1990 (slice bathing solution)
        cao0_ca_ion = 2.4
        ko0_k_ion = 6.24
        nao0_na_ion = 150
        
        # parameters (mho/cm2)
       	my_gna_na = 0.007  # Fit according to AP height (0.02)
       	my_gna_nap = 1.0e-5  # was 1.8e-5 for in vitro testing
       	my_gk_kdrf = 0.004   # Kv3.1 (0.0040)
       	my_gk_iA = 0.001     # Kv4
       	my_gk_skca = 0.0001
       	my_gcaL = 0.002
       	my_gcaN = 0.012
       	
       	# create self.nSoma soma compartments
       	self.soma = [h.Section() for ii in range(self.nSoma)]
 
       	
       	for ii in range(0,int(self.nSoma)):
            self.soma[ii].nseg = 1
            self.soma[ii].Ra = 174
            self.soma[ii].cm = 1
            # Passive Membrane properties
            # self.soma[ii].insert('MyIons') # Dummy mechanism to set up ion concentrations for Na and K
            # Autonomous fast firing membrane properties
            
            self.soma[ii].insert('NaP') # Model for NA+ Leak Current 
            self.soma[ii].insert('Na') # Model for NA+ Current 
            self.soma[ii].insert('KDRf') # Model for K+ Delayed Rectifier Current (fast-deactivating)
            self.soma[ii].insert('iA') # Model for A-current
            self.soma[ii].insert('sKCa')                
            self.soma[ii].insert('cacum') # Model for CA2+ Accumulation
            self.soma[ii].insert('HVA') # Model for HVA Ca2+ Current    
            self.soma[ii].insert('extracellular')            
            h.cai0_cacum = cai0_ca_ion
            for seg in self.soma[ii]:
                seg.gna_Na = my_gna_na
                seg.gnap_NaP = my_gna_nap 
                seg.gk_KDRf = my_gk_kdrf
                seg.gkbar_iA = my_gk_iA
                seg.gk_sKCa = my_gk_skca
                seg.HVA.gcaL = my_gcaL
                seg.HVA.gcaN = my_gcaN
                seg.xraxial[0] = 1e09
                seg.xg[0] = 1e09
                seg.xc[0] = 0


                     
       	# create self.nDend dendrite compartments
       	self.dend = [h.Section() for ii in range(self.nDend)]
       	for ii in range(0,int(self.nDend)):
            self.dend[ii].nseg = 1
            self.dend[ii].Ra = 174
            self.dend[ii].cm = 1
            # Passive Membrane properties
            # self.dend[ii].insert('myions') # Dummy mechanism to set up ion concentrations for Na and K
            # Autonomous fast firing membrane properties
            self.dend[ii].insert('Na') # Model for NA+ Current
            self.dend[ii].insert('NaP') # Model for NA+ Leak Current  
            self.dend[ii].insert('KDRf') # Model for K+ Delayed Rectifier Current (fast-deactivating)
            self.dend[ii].insert('iA') # Model for A-current
            self.dend[ii].insert('sKCa')                
            self.dend[ii].insert('cacum') # Model for CA2+ Accumulation
            self.dend[ii].insert('HVA') # Model for HVA Ca2+ Current    
            self.dend[ii].insert('extracellular')
            h.cai0_cacum = cai0_ca_ion            
            for seg in self.dend[ii]:
                seg.gna_Na = my_gna_na
                seg.gnap_NaP = my_gna_nap
                seg.gk_KDRf = my_gk_kdrf
                seg.gkbar_iA = my_gk_iA
                seg.gk_sKCa = my_gk_skca
                seg.HVA.gcaL = my_gcaL
                seg.HVA.gcaN = my_gcaN
                seg.xraxial[0] = 1e9
                seg.xg[0] = 1e9
                seg.xc[0] = 0
                
                
                
       	# create self.nInitseg initial segment compartment
       	self.initseg = [h.Section() for ii in range(self.nInitseg)]
       	
       	for ii in range(0,int(self.nInitseg)):
            self.initseg[ii].nseg = 1
            self.initseg[ii].Ra = 174
            self.initseg[ii].cm = 1
            # Passive Membrane properties
            # self.initseg[ii].insert('myions') # Dummy mechanism to set up ion concentrations for Na and K
            # Autonomous fast firing membrane properties
            self.initseg[ii].insert('Na') # Model for NA+ Current
            self.initseg[ii].insert('NaP') # Model for NA+ Leak Current  
            self.initseg[ii].insert('KDRf') # Model for K+ Delayed Rectifier Current (fast-deactivating)
            self.initseg[ii].insert('iA') # Model for A-current
            self.initseg[ii].insert('sKCa')                
            self.initseg[ii].insert('cacum') # Model for CA2+ Accumulation
            self.initseg[ii].insert('HVA') # Model for HVA Ca2+ Current    
            self.initseg[ii].insert('extracellular')
            h.cai0_cacum = cai0_ca_ion
            for seg in self.initseg[ii]:
                seg.gna_Na = my_gna_na
                seg.gnap_NaP = my_gna_nap 
                seg.gk_KDRf = my_gk_kdrf
                seg.gkbar_iA = my_gk_iA
                seg.gk_sKCa = my_gk_skca
                seg.HVA.gcaL = my_gcaL
                seg.HVA.gcaN = my_gcaN
                seg.xraxial[0] = 1e09
                seg.xg[0] = 1e09
                seg.xc[0] = 0


        # create self.nNodes node compartments
        self.node = [h.Section() for ii in range(self.nNodes)] 
        
        # create all nodes that are not last node and insert mechanisms
        for ii in range(0,int(self.nNodes)-1):
            self.node[ii].nseg = 1
            self.node[ii].Ra = rhoa/10000
            self.node[ii].cm = 2
            self.node[ii].insert('axnode70')
            self.node[ii].insert('extracellular')
            for seg in self.node[ii]:
                seg.xraxial[0] = Rpn0 # MOhms/cm
                seg.xg[0] = 1e10
                seg.xc[0] = 0
        
        # create last node (passive) and insert mechanisms
        self.node[self.nNodes-1].nseg = 1 # number of segments
        self.node[self.nNodes-1].Ra = rhoa/10000 # axial resistivity (Ohm-cm)
        self.node[self.nNodes-1].cm = 2
        self.node[self.nNodes-1].insert('pas') # insert passive current
        self.node[self.nNodes-1].insert('extracellular') # insert extracellular mechanism
        for seg in self.node[self.nNodes-1]:
            seg.pas.g = 0.0001 # S/cm^2
            seg.pas.e = self.v_init # mV
            seg.xraxial[0] = Rpn0 # MOhms/cm
            seg.xg[0] = 1e10
            seg.xc[0] = 0
        
        
        # create self.nMysa mysa compartments and insert mechanisms
        self.MYSA = [h.Section() for ii in range(0,self.nMysa)] 
        
        for ii in range(0,int(self.nMysa)):
            self.MYSA[ii].nseg = 1
            self.MYSA[ii].Ra = rhoa/10000
            self.MYSA[ii].cm = 2
            self.MYSA[ii].insert('pas')
            self.MYSA[ii].insert('extracellular') # insert extracellular mechanism
            for seg in self.MYSA[ii]:
                seg.pas.g = 0.0001 # S/cm^2
                seg.pas.e = self.v_init # mV
                seg.xraxial[0] = Rpn1 # MOhms/cm
                seg.xg[0] = mygm/(nl*2)
                seg.xc[0] = mycm/(nl*2)	
        
        
        # create self.nFlut flut compartments and insert mechanisms
        self.FLUT = [h.Section() for ii in range(0,self.nFlut)]  
        
        for ii in range(0,int(self.nFlut)):
            self.FLUT[ii].nseg = 1
            self.FLUT[ii].Ra = rhoa/10000
            self.FLUT[ii].cm = 2
            self.FLUT[ii].insert('parak70')
            self.FLUT[ii].insert('pas') # insert passive current
            self.FLUT[ii].insert('extracellular') # insert extracellular mechanism
            for seg in self.FLUT[ii]:
                seg.pas.g = 0.0001 # S/cm^2
                seg.pas.e = self.v_init # mV
                seg.xraxial[0] = Rpn2 # MOhms/cm
                seg.xg[0] = mygm/(nl*2)
                seg.xc[0] = mycm/(nl*2)
        
        
        # create self.nStin stin compartments and insert mechanisms
        self.STIN = [h.Section() for ii in range(0,self.nStin)]  
        for ii in range(0,int(self.nStin)):
            self.STIN[ii].nseg = 1
            self.STIN[ii].Ra = rhoa/10000
            self.STIN[ii].cm = 2
            self.STIN[ii].insert('pas') # insert passive current
            self.STIN[ii].insert('extracellular') # insert extracellular mechanism
            for seg in self.STIN[ii]:
                seg.pas.g = 0.0001 # S/cm^2
                seg.pas.e = self.v_init # mV
                seg.xraxial[0] = Rpx # MOhms/cm
                seg.xg[0] = mygm/(nl*2)
                seg.xc[0] = mycm/(nl*2)
        
        #---------------------------------------------------------------------------
        # Connect the sections in order to make a single multi-compartment cell with soma, dendrites, initial segment, and axons
        #   axon order: mysa, flut, stin, stin, stin, flut, mysa 
        #---------------------------------------------------------------------------
        # connect segments
        # hoc: connect child(cx),parent(px)
        # python: child.connect(parent,px,cx)
        
        # Somas
        for ii in range(self.nSoma-1):
            self.soma[ii+1].connect(self.soma[ii],(1),(0))
            
        self.dend[0].connect(self.soma[2],(1),(0))
        
        #Dendrites
        for ii in range(0,35):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
        
        self.dend[36].connect(self.dend[18](1),(0))
        
        for ii in range(36,51):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[52].connect(self.dend[11](1),(0))    
        
        for ii in range(52,70):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[71].connect(self.dend[61](1),(0))   
        
        for ii in range(71,88):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[89].connect(self.dend[8](1),(0))   
        
        for ii in range(89,117):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[118].connect(self.dend[98](1),(0))   
        
        for ii in range(118,131):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[132].connect(self.dend[8](1),(0))  
        
        for ii in range(132,207):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[208].connect(self.dend[171](1),(0))   
            
        for ii in range(208,234):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[235].connect(self.dend[8](1),(0))   
        
        for ii in range(235,290):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.dend[291].connect(self.dend[259](1),(0))   
        
        for ii in range(291,296):
            self.dend[ii+1].connect(self.dend[ii],(1),(0))
            
        self.initseg[0].connect(self.soma[12](1),(0))
        self.node[0].connect(self.initseg[0](1),(0)) 
        
        #Axons
        #Thalamic axon
        for ii in range(0,48):
            self.MYSA[2*ii].connect(self.node[ii],(1),(0))
            self.FLUT[2*ii].connect(self.MYSA[2*ii],(1),(0))
            
            x = self.nStinPerStretch
            self.STIN[x*ii].connect(self.FLUT[2*ii],(1),(0))
            for gg in range(0,self.nStinPerStretch-1):
                self.STIN[x*ii+gg+1].connect(self.STIN[x*ii+gg],(1),(0))
            
            self.FLUT[2*ii+1].connect(self.STIN[x*ii+x-1],(1),(0))
            self.MYSA[2*ii+1].connect(self.FLUT[2*ii+1],(1),(0))
            self.node[ii+1].connect(self.MYSA[2*ii+1],(1),(0))
        
        self.MYSA[96].connect(self.node[48],(1),(0))
        self.FLUT[96].connect(self.MYSA[96],(1),(0))
        self.STIN[144].connect(self.FLUT[96],(1),(0))
        
        #SNc Axon branch
        self.STIN[145].connect(self.STIN[144],(1),(0))
        self.FLUT[97].connect(self.STIN[145],(1),(0))
        self.MYSA[97].connect(self.FLUT[97],(1),(0))
        self.node[49].connect(self.MYSA[97],(1),(0))
        
        for ii in range(49,78):
            self.MYSA[2*ii].connect(self.node[ii],(1),(0))
            self.FLUT[2*ii].connect(self.MYSA[2*ii],(1),(0))
            
            x = self.nStinPerStretch
            self.STIN[x*ii-1].connect(self.FLUT[2*ii],(1),(0))
            self.STIN[x*ii].connect(self.STIN[x*ii-1],(1),(0))
            self.STIN[x*ii+1].connect(self.STIN[x*ii],(1),(0))
            self.FLUT[2*ii+1].connect(self.STIN[x*ii+1],(1),(0))
            self.MYSA[2*ii+1].connect(self.FLUT[2*ii+1],(1),(0))
            self.node[ii+1].connect(self.MYSA[2*ii+1],(1),(0))
            
        self.MYSA[156].connect(self.node[78],(1),(0))
        self.FLUT[156].connect(self.MYSA[156],(1),(0))
        self.STIN[233].connect(self.FLUT[156],(1),(0))
        
        #Snc Axon branch
        self.FLUT[157].connect(self.FLUT[3],(1),(0))
        self.MYSA[157].connect(self.FLUT[157],(1),(0))
        self.node[79].connect(self.MYSA[157],(1),(0))
        
        for ii in range(79,109):
            self.MYSA[2*ii].connect(self.node[ii],(1),(0))
            self.FLUT[2*ii].connect(self.MYSA[2*ii],(1),(0))
            
            x = self.nStinPerStretch
            self.STIN[x*ii-3].connect(self.FLUT[2*ii],(1),(0))
            self.STIN[x*ii-2].connect(self.STIN[x*ii-3],(1),(0))
            self.STIN[x*ii-1].connect(self.STIN[x*ii-2],(1),(0))
            self.FLUT[2*ii+1].connect(self.STIN[x*ii-1],(1),(0))
            self.MYSA[2*ii+1].connect(self.FLUT[2*ii+1],(1),(0))
            self.node[ii+1].connect(self.MYSA[2*ii+1],(1),(0))

        self.MYSA[218].connect(self.node[109],(1),(0))
        self.FLUT[218].connect(self.MYSA[218],(1),(0))
        self.STIN[324].connect(self.FLUT[218],(1),(0))
        
        #Descending axon branch
        self.STIN[325].connect(self.STIN[4],(1),(0))
        self.STIN[326].connect(self.STIN[325],(1),(0))
        self.FLUT[219].connect(self.STIN[326],(1),(0))
        self.MYSA[219].connect(self.FLUT[219],(1),(0))
        self.node[110].connect(self.MYSA[219],(1),(0))
        
        for ii in range(110,139):
            self.MYSA[2*ii].connect(self.node[ii],(1),(0))
            self.FLUT[2*ii].connect(self.MYSA[2*ii],(1),(0))
            
            x = self.nStinPerStretch
            self.STIN[x*ii-3].connect(self.FLUT[2*ii],(1),(0))
            self.STIN[x*ii-2].connect(self.STIN[x*ii-3],(1),(0))
            self.STIN[x*ii-1].connect(self.STIN[x*ii-2],(1),(0))
            self.FLUT[2*ii+1].connect(self.STIN[x*ii-1],(1),(0))
            self.MYSA[2*ii+1].connect(self.FLUT[2*ii+1],(1),(0))
            self.node[ii+1].connect(self.MYSA[2*ii+1],(1),(0))
        
        self.MYSA[278].connect(self.node[139],(1),(0))
        self.FLUT[278].connect(self.MYSA[278],(1),(0))
        self.STIN[414].connect(self.FLUT[278],(1),(0))
        
        #---------------------------------------------------------------------------
        # Create the actual geometry by placing each compartment at its correct 3-D 
        # coordinate 
        #---------------------------------------------------------------------------
        
        for ii in range(self.nSoma):
            h.pt3dadd(soma[ii,0],soma[ii,1],soma[ii,2],10,sec=self.soma[ii])
            h.pt3dadd(soma[ii,3],soma[ii,4],soma[ii,5],10,sec=self.soma[ii])
        for ii in range(self.nDend):
            h.pt3dadd(dend[ii,0],dend[ii,1],dend[ii,2],2,sec=self.dend[ii])
            h.pt3dadd(dend[ii,3],dend[ii,4],dend[ii,5],2,sec=self.dend[ii])
        if self.nInitseg ==1:
            h.pt3dadd(initseg[0],initseg[1],initseg[2],10,sec=self.initseg[0])
            h.pt3dadd(initseg[3],initseg[4],initseg[5],10,sec=self.initseg[0])
        else:
            for ii in range(self.nInitseg):
                h.pt3dadd(initseg[ii,0],initseg[ii,1],initseg[ii,2],10,sec=self.initseg[ii])
                h.pt3dadd(initseg[ii,3],initseg[ii,4],initseg[ii,5],10,sec=self.initseg[ii])
        for ii in range(self.nNodes):
            h.pt3dadd(node[ii,0],node[ii,1],node[ii,2],nodeD,sec=self.node[ii])
            h.pt3dadd(node[ii,3],node[ii,4],node[ii,5],nodeD,sec=self.node[ii])
        for ii in range(self.nMysa):
            h.pt3dadd(mysa[ii,0],mysa[ii,1],mysa[ii,2],paraD1,sec=self.MYSA[ii])
            h.pt3dadd(mysa[ii,3],mysa[ii,4],mysa[ii,5],paraD1,sec=self.MYSA[ii])
        for ii in range(self.nFlut):
            h.pt3dadd(flut[ii,0],flut[ii,1],flut[ii,2],paraD2,sec=self.FLUT[ii])
            h.pt3dadd(flut[ii,3],flut[ii,4],flut[ii,5],paraD2,sec=self.FLUT[ii])
        for ii in range(self.nStin):
            h.pt3dadd(stin[ii,0],stin[ii,1],stin[ii,2],axonD,sec=self.STIN[ii])
            h.pt3dadd(stin[ii,3],stin[ii,4],stin[ii,5],axonD,sec=self.STIN[ii])
        
        
        # Create a list of all the sections
        self.allseclist = self.create_seclist()
        
        # Create action potential (AP) counters and membrane voltage recorders
        self.insert_APcounters()
       


    #---------------------------------------------------------------------------
    # Create a list of all the sections
    #---------------------------------------------------------------------------
    def create_seclist(self):
        allseclist = h.SectionList()
        for ii in range(self.nSoma):
            allseclist.append(sec=self.soma[ii])
        for ii in range(self.nDend):
            allseclist.append(sec=self.dend[ii])         
        if self.nInitseg ==1:
            allseclist.append(sec=self.initseg[0])
        else:
            for ii in range(self.nInitseg):
                allseclist.append(sec=self.initseg[ii])
        for ii in range(self.nNodes-1):
            allseclist.append(sec=self.node[ii])
        allseclist.append(sec=self.node[self.nNodes-1])
        for ii in range(self.nMysa):
            allseclist.append(sec=self.MYSA[ii])
        for ii in range(self.nFlut):
            allseclist.append(sec=self.FLUT[ii])
        for ii in range(self.nStin):
            allseclist.append(sec=self.STIN[ii]) 
        return allseclist
            
    #---------------------------------------------------------------------------
    # Create action potential (AP) counters and membrane voltage recorders
    #---------------------------------------------------------------------------
    def insert_APcounters(self):
        self.recNode = 47
        # Create a single action potential counter with a "thresh" setting in order 
        # to capture action potentials and determine the threshold for the axon
        self.apc = h.APCount(0.5,sec=self.node[self.recNode]) 
        self.apc_times = h.Vector()     # empty vector for action potential times
        self.apc.thresh = -20           # threshold that must be crossed to count as an action potential (mV)
        self.apc.record(self.apc_times) # set it to record

        # Set up AP counters in every node and branches (if axon has any)
        self.nAPcounters = self.nNodes
        
        # records AP time at every node
        self.AP_timerec = []
        for ii in range(self.nAPcounters):
            self.AP_timerec.append(h.Vector())
    
        # Setup AP counter
        self.AP_counters = []
        for ii in range(self.nNodes):
            self.AP_counters.append(h.APCount(0.5,sec=self.node[ii]))    # put AP counter in every node
            self.AP_counters[ii].record(self.AP_timerec[ii])                  # records times of APs
        
        # Setup membrane voltage recorder at each time step
        self.membrane_v = h.Vector()
        self.membrane_v_time = h.Vector()
        self.membrane_v.record(self.node[self.recNode](0.5)._ref_v)
        self.membrane_v_time.record(h._ref_t)
        
    
class Axon(NRTL):
    """
    This Axon class inherits the NRTL metaclass methods.
       
    Methods:     
        construct_geom - build axon geometry (also calls create_seclist and insert_APcounters
        create_seclist - creates a list of all the sections
        insert_APcounters - insert action potential counters and membrane voltage recorders
        
    Parent (NRTL) Methods:     
        simulation initialization - passes simulation relevant attributes from a StimulusTrain object
        insert_vext - calculated and plays the extracellular voltages into the multicompartment neuron
        integrate - solve the multicompartment neuron equivalent circuit model at each time point using NERUON 
        response_to_stimulus - determine if a neuron was Activated
        trial - insert_vext, integrate, and response_to_stimulus in that order
        threshold - determine the neuron threshold for Activation by running multiple trials
        
    NEURON Globals:
        h.dt:           (double) integration time-step for neuron (msec)
        h.celsius:      (double) temperature in Celsius
        
    Attributes:
        v_init:          (double) cell initial voltage
        nStinPerStretch: (int) number of STIN in a row
        nNodes:          (int) number of comparmtents of type "node"
        nMysa:           (int) number of comparmtents of type "mysa"
        nFlut:           (int) number of comparmtents of type "flut"
        nStin:           (int) number of comparmtents of type "stin"
        node:            (NEURON Section object) one for each "node" compartment
        MYSA:            (NEURON Section object) one for each "mysa" compartment
        FLUT:            (NEURON Section object) one for each "flut" compartment
        STIN:            (NEURON Section object) one for each "stin" compartment
        allseclist       (NERUON list) of sections
        apc:             (NEURON APCount object) for counter APs
        apc_times:       (NEURON Vector object) for storing times of APs
        nAPcounters:     (int) number of AP counters, equal to nNodes          
        AP_counters:     (list of NEURON APCount objects) one on each node
        AP_timerec:      (NEURON Vector object) for storing times of APs
        membrane_v:      (NEURON Vector object) for recording the membrane voltage on the last node
        membrane_v_time: (NEURON Vector object) for storing the timepoints related to the membrane voltage on the last node
    """
    
    def construct_geom(self,axon_data,fiberD=2):        
        """
        Build a multicompartment axon, generate a section list, and insert AP 
        counters 
        
        >>> axon_obj.construct_geom(axon_data,fiberD)
        
        Inputs:
            fiberD:    (double) fiber diameter (um): must be 1, 2 (default), 5.7, 7.3, or 8.7
            axon_data: (structure) imported from a .mat and contains all the compartment information
        """   
        # Set axon parameters and build axons
        h.celsius = 36  # temperature in Celsius
        self.v_init = -70  # axon resting potential
      
        # get axon compartment data
        nodes = axon_data.nodes
        mysa = axon_data.mysa
        flut = axon_data.flut
        stin = axon_data.stin
        self.nStinPerStretch = axon_data.nStinPerStretch
        
        # calculate number of each compartment type
        self.nNodes = int(nodes.shape[0]) # number of nodes
        self.nMysa = int(mysa.shape[0]) # number of mysa
        self.nFlut = int(flut.shape[0]) # number of flut
        self.nStin = int(stin.shape[0]) # number of internodes
        total = self.nNodes+self.nMysa+self.nFlut+self.nStin # total number of compartments
    
        # Axon electrical parameters
        rhoa=0.7e6 # Ohm-um
        mycm=0.1   # uF/cm2/lamella membrane
        mygm=0.001 # S/cm2/lamella membrane
    
    
        #-----------------------------------------------------------------------
        # Set morphological parameters based on fiber diameter
        #-----------------------------------------------------------------------
        space_p1 = 0.002  
        space_p2 = 0.004
        space_i = 0.004
        if fiberD==1: axonD=0.8; nodeD=0.7; paraD1=0.7; paraD2=0.8; nl=20
        elif fiberD==2: axonD=1.6; nodeD=1.4; paraD1=1.4; paraD2=1.6; nl=30
        elif fiberD==5.7: axonD=3.4; nodeD=1.9; paraD1=1.9; paraD2=3.4; nl=80
        elif fiberD==7.3: axonD=4.6; nodeD=2.4; paraD1=2.4; paraD2=4.6; nl=100
        elif fiberD==8.7: axonD=5.8; nodeD=2.8; paraD1=2.8; paraD2=5.8; nl=110
        else: raise Exception("Unknown fiber diameter")
    
        Rpn0=(rhoa*.01)/(pi*((((nodeD/2)+space_p1)**2)-((nodeD/2)**2)))
        Rpn1=(rhoa*.01)/(pi*((((paraD1/2)+space_p1)**2)-((paraD1/2)**2)))
        Rpn2=(rhoa*.01)/(pi*((((paraD2/2)+space_p2)**2)-((paraD2/2)**2)))
        Rpx=(rhoa*.01)/(pi*((((axonD/2)+space_i)**2)-((axonD/2)**2)))
    
    
        #-----------------------------------------------------------------------
        # Create "sections" for each compartment and assign mechanisms
        #-----------------------------------------------------------------------
        # create self.nNodes node compartments
        self.node = [h.Section() for ii in range(self.nNodes)] 
        
        # create first node (passive) and insert mechanisms
        self.node[0].nseg = 1 # number of segments
        self.node[0].Ra = rhoa/10000 # axial resistivity (Ohm-cm)
        self.node[0].cm = 2
        self.node[0].insert('pas') # insert passive current
        self.node[0].insert('extracellular') # insert extracellular mechanism
        for seg in self.node[0]:
            seg.pas.g = 0.0001 # S/cm^2
            seg.pas.e = self.v_init # mV
            seg.xraxial[0] = Rpn0 # MOhms/cm
            seg.xg[0] = 1e10
            seg.xc[0] = 0
    
        # create all nodes that are not first or last node and insert mechanisms
        for ii in range(1,self.nNodes-1):
            self.node[ii].nseg = 1
            self.node[ii].Ra = rhoa/10000
            self.node[ii].cm = 2
            self.node[ii].insert('axnode70')
            self.node[ii].insert('extracellular')
            for seg in self.node[ii]:
                seg.xraxial[0] = Rpn0 # MOhms/cm
                seg.xg[0] = 1e10
                seg.xc[0] = 0
    
        # create last node (passive) and insert mechanisms
        self.node[self.nNodes-1].nseg = 1 # number of segments
        self.node[self.nNodes-1].Ra = rhoa/10000 # axial resistivity (Ohm-cm)
        self.node[self.nNodes-1].cm = 2
        self.node[self.nNodes-1].insert('pas') # insert passive current
        self.node[self.nNodes-1].insert('extracellular') # insert extracellular mechanism
        for seg in self.node[self.nNodes-1]:
            seg.pas.g = 0.0001 # S/cm^2
            seg.pas.e = self.v_init # mV
            seg.xraxial[0] = Rpn0 # MOhms/cm
            seg.xg[0] = 1e10
            seg.xc[0] = 0
   	
   	
        # create self.nMysa mysa compartments and insert mechanisms
        self.MYSA = [h.Section() for ii in range(0,self.nMysa)] 
        for ii in range(0,int(self.nMysa)):
            self.MYSA[ii].nseg = 1
            self.MYSA[ii].Ra = rhoa/10000
            self.MYSA[ii].cm = 2
            self.MYSA[ii].insert('pas')
            self.MYSA[ii].insert('extracellular') # insert extracellular mechanism
            for seg in self.MYSA[ii]:
                seg.pas.g = 0.0001 # S/cm^2
                seg.pas.e = self.v_init # mV
                seg.xraxial[0] = Rpn1 # MOhms/cm
                seg.xg[0] = mygm/(nl*2)
                seg.xc[0] = mycm/(nl*2)	
    
    
        # create self.nFlut flut compartments and insert mechanisms
        self.FLUT = [h.Section() for ii in range(0,self.nFlut)]  
        for ii in range(0,int(self.nFlut)):
            self.FLUT[ii].nseg = 1
            self.FLUT[ii].Ra = rhoa/10000
            self.FLUT[ii].cm = 2
            self.FLUT[ii].insert('parak70')
            self.FLUT[ii].insert('pas') # insert passive current
            self.FLUT[ii].insert('extracellular') # insert extracellular mechanism
            for seg in self.FLUT[ii]:
                seg.pas.g = 0.0001 # S/cm^2
                seg.pas.e = self.v_init # mV
                seg.xraxial[0] = Rpn2 # MOhms/cm
                seg.xg[0] = mygm/(nl*2)
                seg.xc[0] = mycm/(nl*2)
    
   	    
        # create self.nStin stin compartments and insert mechanisms
        self.STIN = [h.Section() for ii in range(0,self.nStin)]  
        for ii in range(0,int(self.nStin)):
            self.STIN[ii].nseg = 1
            self.STIN[ii].Ra = rhoa/10000
            self.STIN[ii].cm = 2
            self.STIN[ii].insert('pas') # insert passive current
            self.STIN[ii].insert('extracellular') # insert extracellular mechanism
            for seg in self.STIN[ii]:
                seg.pas.g = 0.0001 # S/cm^2
                seg.pas.e = self.v_init # mV
                seg.xraxial[0] = Rpx # MOhms/cm
                seg.xg[0] = mygm/(nl*2)
                seg.xc[0] = mycm/(nl*2)
    
    
        #---------------------------------------------------------------------------
        # Connect the sections in order to make a single multi-compartment axon
        #     order: mysa, flut, stin, stin, stin, flut, mysa 
        #---------------------------------------------------------------------------
        # connect segments
        # hoc: connect child(cx),parent(px)
        # python: child.connect(parent,px,cx)
        for ii in range(self.nNodes-1):
            self.MYSA[2*ii].connect(self.node[ii],(1),(0))
            self.FLUT[2*ii].connect(self.MYSA[2*ii],(1),(0))
  		
            x = self.nStinPerStretch
            self.STIN[x*ii].connect(self.FLUT[2*ii],(1),(0))
            for gg in range(0,self.nStinPerStretch-1):
                self.STIN[x*ii+gg+1].connect(self.STIN[x*ii+gg],(1),(0))
    
            self.FLUT[2*ii+1].connect(self.STIN[x*ii+x-1],(1),(0))
            self.MYSA[2*ii+1].connect(self.FLUT[2*ii+1],(1),(0))
            self.node[ii+1].connect(self.MYSA[2*ii+1],(1),(0))
           	
            # depending on the number of points left over connect according to the 
            # predetermined organization of the axon design: 
            #mysa, flut, stin, stin, stin, flut, mysa 
        
        remainingpts = total - ((self.nStinPerStretch+5)*(self.nNodes-1) + 1)
        if remainingpts and self.nStinPerStretch==3: # if there are left-over points
            if remainingpts==1:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
            
            elif remainingpts==2:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
            
            elif remainingpts==3:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-3].connect(self.FLUT[self.nFlut-1],(1),(0))
            
            elif remainingpts==4:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-3].connect(self.FLUT[self.nFlut-1],(1),(0))
                self.STIN[self.nStin-2].connect(self.STIN[self.nStin-3],(1),(0))
            
            elif remainingpts==5:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-3].connect(self.FLUT[self.nFlut-1],(1),(0))
                self.STIN[self.nStin-2].connect(self.STIN[self.nStin-3],(1),(0))
                self.STIN[self.nStin-1].connect(self.STIN[self.nStin-2],(1),(0))
            
            elif remainingpts==6:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-2].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-3].connect(self.FLUT[self.nFlut-2],(1),(0))
                self.STIN[self.nStin-2].connect(self.STIN[self.nStin-3],(1),(0))
                self.STIN[self.nStin-1].connect(self.STIN[self.nStin-2],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.STIN[self.nStin-1],(1),(0))
            
            elif remainingpts==7:
                self.MYSA[self.nMysa-2].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-2].connect(self.MYSA[self.nMysa-2],(1),(0))
                self.STIN[self.nStin-3].connect(self.FLUT[self.nFlut-2],(1),(0))
                self.STIN[self.nStin-2].connect(self.STIN[self.nStin-3],(1),(0))
                self.STIN[self.nStin-1].connect(self.STIN[self.nStin-2],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.STIN[self.nStin-1],(1),(0))
                self.MYSA[self.nMysa-1].connect(self.FLUT[self.nFlut-1],(1),(0))
            else:
                raise Exception("More than 7 points remaining")
        elif remainingpts and self.nStinPerStretch==6: 
            if remainingpts==1:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
            
            elif remainingpts==2:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
            
            elif remainingpts==3:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-6].connect(self.FLUT[self.nFlut-1],(1),(0))
            
            elif remainingpts==4:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-6].connect(self.FLUT[self.nFlut-1],(1),(0))
                self.STIN[self.nStin-5].connect(self.STIN[self.nStin-6],(1),(0))
            
            elif remainingpts==5:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-6].connect(self.FLUT[self.nFlut-1],(1),(0))
                self.STIN[self.nStin-5].connect(self.STIN[self.nStin-6],(1),(0))
                self.STIN[self.nStin-4].connect(self.STIN[self.nStin-5],(1),(0))
            
            elif remainingpts==6:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-6].connect(self.FLUT[self.nFlut-1],(1),(0))
                self.STIN[self.nStin-5].connect(self.STIN[self.nStin-6],(1),(0))
                self.STIN[self.nStin-4].connect(self.STIN[self.nStin-5],(1),(0))
                self.STIN[self.nStin-3].connect(self.STIN[self.nStin-4],(1),(0))
                                
            elif remainingpts==7:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-6].connect(self.FLUT[self.nFlut-1],(1),(0))
                self.STIN[self.nStin-5].connect(self.STIN[self.nStin-6],(1),(0))
                self.STIN[self.nStin-4].connect(self.STIN[self.nStin-5],(1),(0))
                self.STIN[self.nStin-3].connect(self.STIN[self.nStin-4],(1),(0))
                self.STIN[self.nStin-2].connect(self.STIN[self.nStin-3],(1),(0))
            
            elif remainingpts==8:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-6].connect(self.FLUT[self.nFlut-1],(1),(0))
                self.STIN[self.nStin-5].connect(self.STIN[self.nStin-6],(1),(0))
                self.STIN[self.nStin-4].connect(self.STIN[self.nStin-5],(1),(0))
                self.STIN[self.nStin-3].connect(self.STIN[self.nStin-4],(1),(0))
                self.STIN[self.nStin-2].connect(self.STIN[self.nStin-3],(1),(0))
                self.STIN[self.nStin-1].connect(self.STIN[self.nStin-2],(1),(0))
            
            elif remainingpts==9:
                self.MYSA[self.nMysa-1].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-2].connect(self.MYSA[self.nMysa-1],(1),(0))
                self.STIN[self.nStin-6].connect(self.FLUT[self.nFlut-2],(1),(0))
                self.STIN[self.nStin-5].connect(self.STIN[self.nStin-6],(1),(0))
                self.STIN[self.nStin-4].connect(self.STIN[self.nStin-5],(1),(0))
                self.STIN[self.nStin-3].connect(self.STIN[self.nStin-4],(1),(0))
                self.STIN[self.nStin-2].connect(self.STIN[self.nStin-3],(1),(0))
                self.STIN[self.nStin-1].connect(self.STIN[self.nStin-2],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.STIN[self.nStin-1],(1),(0))
                
            
            elif remainingpts==10:
                self.MYSA[self.nMysa-2].connect(self.node[self.nNodes-1],(1),(0))
                self.FLUT[self.nFlut-2].connect(self.MYSA[self.nMysa-2],(1),(0))
                self.STIN[self.nStin-6].connect(self.FLUT[self.nFlut-2],(1),(0))
                self.STIN[self.nStin-5].connect(self.STIN[self.nStin-6],(1),(0))
                self.STIN[self.nStin-4].connect(self.STIN[self.nStin-5],(1),(0))
                self.STIN[self.nStin-3].connect(self.STIN[self.nStin-4],(1),(0))
                self.STIN[self.nStin-2].connect(self.STIN[self.nStin-3],(1),(0))
                self.STIN[self.nStin-1].connect(self.STIN[self.nStin-2],(1),(0))
                self.FLUT[self.nFlut-1].connect(self.STIN[self.nStin-1],(1),(0))
                self.MYSA[self.nMysa-1].connect(self.FLUT[self.nFlut-1],(1),(0))
            else:
                raise Exception("More than 10 points remaining")
        elif remainingpts != 0: 
            raise Exception("Number of Stins per Segment not supported") 
        #---------------------------------------------------------------------------
        # Create the actual geometry by placing each compartment at its correct 3-D 
        # coordinate 
        #---------------------------------------------------------------------------
        for ii in range(self.nNodes):
            h.pt3dadd(nodes[ii,0],nodes[ii,1],nodes[ii,2],nodeD,sec=self.node[ii])
            h.pt3dadd(nodes[ii,3],nodes[ii,4],nodes[ii,5],nodeD,sec=self.node[ii])
        for ii in range(self.nMysa):
            h.pt3dadd(mysa[ii,0],mysa[ii,1],mysa[ii,2],paraD1,sec=self.MYSA[ii])
            h.pt3dadd(mysa[ii,3],mysa[ii,4],mysa[ii,5],paraD1,sec=self.MYSA[ii])
        for ii in range(self.nFlut):
            h.pt3dadd(flut[ii,0],flut[ii,1],flut[ii,2],paraD2,sec=self.FLUT[ii])
            h.pt3dadd(flut[ii,3],flut[ii,4],flut[ii,5],paraD2,sec=self.FLUT[ii])
        for ii in range(self.nStin):
            h.pt3dadd(stin[ii,0],stin[ii,1],stin[ii,2],axonD,sec=self.STIN[ii])
            h.pt3dadd(stin[ii,3],stin[ii,4],stin[ii,5],axonD,sec=self.STIN[ii])
        
        
        # Create a list of all the sections
        self.allseclist = self.create_seclist()
        
        # Create action potential (AP) counters and membrane voltage recorders
        self.insert_APcounters()
        

    #---------------------------------------------------------------------------
    # Create a list of all the sections
    #---------------------------------------------------------------------------
    def create_seclist(self):
        allseclist = h.SectionList()
        for ii in range(self.nNodes):
            allseclist.append(sec=self.node[ii])
        for ii in range(self.nMysa):
            allseclist.append(sec=self.MYSA[ii])
        for ii in range(self.nFlut):
            allseclist.append(sec=self.FLUT[ii])
        for ii in range(self.nStin):
            allseclist.append(sec=self.STIN[ii]) 
        return allseclist
            
    #---------------------------------------------------------------------------
    # Create action potential (AP) counters and membrane voltage recorders
    #---------------------------------------------------------------------------
    def insert_APcounters(self):
        self.recNode = self.nNodes-1 #0
        # Create a single action potential counter with a "thresh" setting in order 
        # to capture action potentials and determine the threshold for the axon
        self.apc = h.APCount(0.5,sec=self.node[self.recNode])  # always place this AP counter on the last node so it does not count APs that blocked due to collision (only really important for axons that fire passively)
        self.apc_times = h.Vector()     # empty vector for action potential times
        self.apc.thresh = -20           # threshold that must be crossed to count as an action potential (mV)
        self.apc.record(self.apc_times) # set it to record

        # Set up AP counters in every node and branches (if axon has any)
        self.nAPcounters = self.nNodes
        
        # records AP time at every node
        self.AP_timerec = []
        for ii in range(self.nAPcounters):
            self.AP_timerec.append(h.Vector())
    
        # Setup AP counter
        self.AP_counters = []
        for ii in range(self.nNodes):
            self.AP_counters.append(h.APCount(0.5,sec=self.node[ii]))    # put AP counter in every node
            self.AP_counters[ii].record(self.AP_timerec[ii])                  # records times of APs
        
        # Setup membrane voltage recorder at each time step
        self.membrane_v = h.Vector()
        self.membrane_v_time = h.Vector()
        self.membrane_v.record(self.node[self.recNode](0.5)._ref_v)
        self.membrane_v_time.record(h._ref_t)

