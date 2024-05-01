% NOTES:
% you may need to route to the dtifit folder of your patient
% make sure file is named: cortex_parcellation-label.nii

user = 'butle821';
disorder = 'Essential_Tremor'; 
patient = 'ET046';

%% Load DWI image and dimensions

dtifit_folder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\' disorder '\' patient '\DTI\dtifit\'];
brain_segmented_by_fast_DWI_space = ['C:\Users\' user '\Box\DBS_Patient_Modeling\' disorder '\' patient '\cortex_parcellation/cortex_parcellation-label.nii'];

% Obtain DWI dimensions from Amira/mricron (open image derived from DWI, e.g. FA).

if strcmp(patient,'ET043')==1, dimensions_of_DWI_images = [168,168,112];
elseif strcmp(patient,'ET044')==1, dimensions_of_DWI_images = [168,168,112]; % in voxels
elseif strcmp(patient,'ET046')==1, dimensions_of_DWI_images = [168,168,112]; % in voxels
elseif strcmp(patient,'ET037')==1, dimensions_of_DWI_images = [168,168,112]; % in voxels 
elseif strcmp(patient,'ET040')==1, dimensions_of_DWI_images = [168,168,112]; % in voxels 
elseif strcmp(patient,'ET025')==1, dimensions_of_DWI_images = [136,136,82]; % in voxels 
elseif strcmp(patient,'ET024')==1, dimensions_of_DWI_images = [136,136,76]; % in voxels
elseif strcmp(patient,'ET023')==1, dimensions_of_DWI_images = [136,136,76]; % in voxels
elseif strcmp(patient,'ET039')==1, dimensions_of_DWI_images = [168,168,112]; % in voxels        
elseif strcmp(patient,'ET036')==1, dimensions_of_DWI_images = [136,136,80]; % in voxels
elseif strcmp(patient,'ET026')==1, dimensions_of_DWI_images = [168,168,112]; % in voxels
elseif strcmp(patient,'ET033')==1, dimensions_of_DWI_images = [136,136, 82]; % in voxels
elseif strcmp(patient,'ET034')==1, dimensions_of_DWI_images = [136,136, 80]; % in voxels
elseif strcmp(patient,'ET035')==1, dimensions_of_DWI_images = [136,136, 78]; % in voxels
end 

dti = 'dtiwls_54dir_S0.nii.gz';
dtifile = ['C:\Users\' user '\Box\DBS_Patient_Modeling\' disorder '\' patient '\DTI\dtifit\' dti];
info = niftiinfo(dtifile);
base_transform_to_DWI_space_on_Slicer = info.Transform.T';
% You may need to use getTransform to find the base transform:
% EX: base_transform_to_DWI_space_on_Amira = reshape([1.5 0 0 0 0 -1.47092 0.293919 0 0 -0.293919 -1.47092 0 -103.569 85.8276 48.9382 1],4,4);

% if strcmp(patient,'ET030') ==1 
%     base_transform_to_DWI_space_on_Slicer = [1.5 0 -0.00051 -98.60710; -0.00014 1.44403 -0.40591 -57.32828; 0.00049 0.40591 1.44403 -86.94191] 
% elseif strcmp(patient,'ET025') ==1 
%     base_transform_to_DWI_space_on_Slicer = [1.49787 0 0.07898 -98.30307; 0.00989 1.48847 -0.18537 -63.53344; -0.07928 0.18564 1.48636 -55.46800] 
% elseif strcmp(patient,'ET024') ==1 
%     base_transform_to_DWI_space_on_Slicer = [1.49498 -0.06527 0.10380 -96.38893; 0.09181 1.43758 -0.41825 -71.69290; -0.08128 0.42320 1.43677 -36.47486] 
% elseif strcmp(patient,'ET023') ==1 
%     base_transform_to_DWI_space_on_Slicer = [1.5 0 0.00001 -101.65289; 0 1.43748 -0.42853 -51.81146; -0.00001 0.42853 1.43748 -53.51049] %MRIcron
% end 

output_folder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\' disorder '\' patient '\cortex_parcellation/'];
is_save_file = 1;

%% Set median frequency & corresponding tissue conductivities

isotropic_electric_properties.freq = 3049; % median freq in [Hz]
% Use 2415 Hz for Abbott pulse delay study, 3049 Hz for current-controlled waveform 

%Calculate white matter & grey matter properties based on the median freq of your waveform using this caluclator: http://niremf.ifac.cnr.it/tissprop/htmlclie/htmlclie.php
%The encapsulation layer should be modeled using white matter properties.
if isotropic_electric_properties.freq == 2415
	isotropic_electric_properties.gray_sigma_iso = 0.10407; % S/m
	isotropic_electric_properties.white_sigma_iso = 0.064443; % S/m 
    isotropic_electric_properties.gray_epsilon_iso = 80409; % S/mf
	isotropic_electric_properties.white_epsilon_iso = 35136; % S/m
elseif isotropic_electric_properties.freq == 3049
    isotropic_electric_properties.gray_sigma_iso = 0.10577; % S/m
    isotropic_electric_properties.white_sigma_iso = 0.065095; % S/m 
    isotropic_electric_properties.gray_epsilon_iso = 65898; % S/mf
    isotropic_electric_properties.white_epsilon_iso = 29790; % S/m
end

%% Generate conductivity tensors

generate_conductivity_tensors_v5(dtifit_folder, brain_segmented_by_fast_DWI_space,...
    dimensions_of_DWI_images, base_transform_to_DWI_space_on_Slicer,output_folder,...
    is_save_file,isotropic_electric_properties)