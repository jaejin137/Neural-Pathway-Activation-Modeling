%function run_simulation_FourierFEM(nhp_name, geom_name, fem_name, electrodes, conductance_model, freq, ETI, ppn, nodes, pmem, walltime)

%% get livelink rolling
% pause(60); % pause for 30 sec to enable comsol livelink with matlab to launch peacefully
% addpath ~/comsol54/multiphysics/mli; % add the path to comsol's mph (livelink) functions
% %addpath /panfs/roc/itascasoft/comsol/5.2a/mli % if using MSI license (although i havent needed this-annie) 
% mphstart; % connect to comsol server
% pause(60); % pause for 30 sec to enable connection to comsol server occur peacefully

% nhp_name = {'FEM_refinement'};
%geom_name={'IC','GPi_Output','PM','SMA','PutGPeGPi','PutGPe','GPeSTNSN','GPeGPiSTNSN'};
user='jaejin';
%geom_name={'GPeSTNSN','GPeGPiSTNSN','PutGPe','PutGPeGPi','GPi_Output'};
geom_name={'GPeSTN','STNMCPath'};
full_ID= 'Barb';
fem_name = 'Barb';
side = 'R';
% electrodes = 1:32;
% conductance_model = {'Schmidt','Wolters','Tuch'};
freq = 3049;    %3049 for Abbott, 4294 for Medtronic
% freq = [4400, 5000000];

% ETI = true;

%ppn = 12;              % processors per node
%nodes = 10;            % number of nodes requested
%pmem = '2580mb';       % mem requested per processor
%walltime = '2:00:00'; % job walltime hh:mm:ss check MSI docs for max

FFEM = true;             % if FFEM use different Neuron code
nFreq = length(freq);
Current_density=1e-3/3.115e-8;


% Folder setup
dataproc_dir = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing';
fem_dir = [dataproc_dir filesep 'COMSOL_Output'];
geom_dir = [dataproc_dir filesep 'Axons' filesep 'Cleaned'];
%load sinusoid_amplitudes_pi_6.mat
%n_amps=size(x1,2);

% import comsol untilities
import com.comsol.model.*
import com.comsol.model.util.*

%modelname = 'Barb_R_STNDBS_rev1_ver3_solved';
modelname = 'Barb_R_STNGPeDBS_rev1_ver3_STN-1-3-0p45_solved';

sim_dir = [dataproc_dir filesep 'Axons' filesep 'Interp_' modelname]; 
mkdir(sim_dir)


for j = 1:length(geom_name) 
    % load axon geometries
    load([geom_dir filesep fem_name '_' side '_' geom_name{j} '.mat'],'axonpop')
    nAxon = size(axonpop,2); % number of cells
    compartment_xyz=[];
    % Create a single large 3xN matrix with all the axon compartments.
    % This greatly increases the speed of interpolation from the FEM
    % solution.
    axon_count=0;
    axon_num=[];
    for iAxon = 1:nAxon % for each axon
        axon_count = axon_count + 1;
        temp = [axonpop(iAxon).midpts(:,1),axonpop(iAxon).midpts(:,2),axonpop(iAxon).midpts(:,3)]'; % compartment coordinates
        compartment_xyz = cat(2,compartment_xyz,temp/10^3); % convert compartment coordinates from um to mm (comsol model units are mm)
        axon_num = cat(2,axon_num,axon_count*ones(1,length(temp)));
    end
    compartment_V = zeros(nFreq,length(compartment_xyz),2);
    
    % modify comsol model and solve
    for i = 1:8
        
        % load comsol model
        tic
        disp('Loading model')
        %model = mphload([fem_dir '\' fem_name '_' side '_elec' num2str(i) '_solved.mph']); % load model
        model = mphload([fem_dir filesep modelname '.mph']); % load model
        disp('Model loading complete')
        toc
        
        %model.physics('ec').feature('ncd1').set('nJ', {num2str(Current_density*x1(1,amp))});
        %model.physics('ec').feature('ncd2').set('nJ', {num2str(Current_density*x1(2,amp))});
        %model.physics('ec').feature('ncd3').set('nJ', {num2str(Current_density*x1(3,amp))});
        
        % create a directory for the comsol log file
        
    %      	tic
    %          disp('Loading model')
    %          model = mphload([sim_dir '/fem_out' i '.mph']); % load model
    %          disp('Model loading complete')
    %         toc
        
        
        % save comsol progress in a log file
        comsol_log = [sim_dir filesep 'fem.log'];
        ModelUtil.showProgress(comsol_log);
        
        % solve model
        %disp('Solving model'); tic
        %model.sol('sol3').runAll;
        %disp('Model solving complete'); toc
        
        % save model
        %disp('Saving model'); tic
        %mphsave(model,[fem_dir '/' fem_name '_' num2str(i) '_fourier_solved.mph'])
        %disp('Model save complete'); toc
        
        % interpolate and save FEM results
        
        disp('Interpolating voltage values...');
    
        
        tic
        for ii = 1%:nFreq
            ii;
            % comsol requres that real and imag parts of the solution be interpolated separately.
            %compartment_V_real = mphinterp(model,'real(V2)','Coord',compartment_xyz,'solnum',ii);
            %compartment_V_imag = mphinterp(model,'imag(V2)','Coord',compartment_xyz,'solnum',ii);
            compartment_V_real = mphinterp(model,'real(V)','Coord',compartment_xyz,'solnum',ii);
            compartment_V_imag = mphinterp(model,'imag(V)','Coord',compartment_xyz,'solnum',ii);
            compartment_V(ii,:,i) = complex(compartment_V_real,compartment_V_imag);
            
        end
        toc
        
        ModelUtil.remove('model'); % remove model
        clearvars model
    end
    disp('Saving interpolation results')
    results_dir = [sim_dir filesep fem_name '_' side '_' geom_name{j} ];
    mkdir(results_dir)
    for iAxon=1:nAxon
        axonpop(iAxon).voltage = compartment_V(:,axon_num==iAxon,:); % append voltage values to axonpop structure
        axon = axonpop(iAxon);
        save([results_dir '/axon' num2str(iAxon) '.mat'], 'axon'); % save workspace
    end
    disp('Finished Saving interpolation results')
end 

% Create and submit pbs job for running Neuron simulations
%nrn_sim_dir = [sim_dir '/' geom_name{iGeom}];
%run_neuron(nrn_sim_dir, ppn, nodes, pmem, walltime, FFEM)



%disp('Saving interpolation results is complete')
% Create and submit pbs job for running Neuron simulations
%run_neuron(sim_dir, ppn, nodes, pmem, walltime, FFEM)


