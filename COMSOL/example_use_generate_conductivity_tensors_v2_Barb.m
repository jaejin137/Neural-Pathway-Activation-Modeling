%% Program to transform FAST input into conductivity tensors for COMSOL
% Edit these first 6 variables to fit your subject and folder structure
% user = 'elecy'
% disorder = 'Globus_Pallidus' 
% subject = 'UD1028_PD104' 
% ID='UD1028'
user = 'jaejin'
disorder = 'PD'
subject = 'Barb'
ID = 'Barb'
%dimensions_of_DWI_images = [168 168 112]; % in voxels
dimensions_of_DWI_images = [212 132 64]
ABT = 1; %0 for Medtronic, 1 for Abbott
%dtifit_folder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\' disorder '\' subject '\DTI\dtifit\'];
dtifit_folder = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Reoriented_Data\diffusion\dti'
%brain_segmented_by_fast_DWI_space = ['C:\Users\' user '\Box\DBS_Patient_Modeling\' disorder '\' subject '\Data Processing\COMSOL_Input\' ID '_FAST_seg_DTIspace_expanded-label.nii.gz'];
%brain_segmented_by_fast_DWI_space = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Moddeling_Data\COMSOL_Input\FAST_Segmentation-CSF_expanded-label.nii.gz'
brain_segmented_by_fast_DWI_space = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\COMSOL Input\ver2\FAST_seg_DTIspace_expanded.nii.gz'
% find this by going to the volume module in slicer and check the image dimension under volume information tab of the FA image
comsol_code_folder = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\COMSOL_Code'
cd (dtifit_folder)
%dti = 'dtiwls_54dir_S0.nii.gz';
dti = 'dti_S0.nii.gz'
info = niftiinfo(dti);
base_transform_to_DWI_space_on_Slicer = info.Transform.T';
% output_folder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\' disorder '\' subject '\Data Processing\COMSOL_Input'];
output_folder = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\COMSOL_Input\ver2'
is_save_file = 1;

%Calculate your white matter and grey matter properties based on the median
%frequency of you waveform using this caluclator: http://niremf.ifac.cnr.it/tissprop/htmlclie/htmlclie.php
%Your encapsulation layer should also be modeled using white matter
%properties 
if ABT
    %parameters for median freq = 3049 (current controlled waveform) 
    isotropic_electric_properties.freq = 3049; %Hz - 3049 Hz median freq for current controlled waveform 

    isotropic_electric_properties.gray_sigma_iso = 0.10577; % S/m
    isotropic_electric_properties.white_sigma_iso = 0.065095; % S/m 

    isotropic_electric_properties.gray_epsilon_iso = 65898; % S/mf
    isotropic_electric_properties.white_epsilon_iso = 29790; % S/m
else
    %parameters for median freq = 4294 (voltage controlled waveform)
    isotropic_electric_properties.freq = 4294;

    isotropic_electric_properties.gray_sigma_iso = 0.10836; % S/m
    isotropic_electric_properties.white_sigma_iso = 0.066184; % S/m 

    isotropic_electric_properties.gray_epsilon_iso = 48617; % S/m
    isotropic_electric_properties.white_epsilon_iso = 23361; % S/m
end

% cd(['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\COMSOL_Code'])
cd(comsol_code_folder)
generate_conductivity_tensors_v5_Barb(dtifit_folder, brain_segmented_by_fast_DWI_space,...
    dimensions_of_DWI_images, base_transform_to_DWI_space_on_Slicer,output_folder,...
    is_save_file,isotropic_electric_properties)

%look at base transform

