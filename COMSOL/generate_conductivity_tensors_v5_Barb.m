% generateConductivityMaps.m
%
% dependency: load_nifti
%
% Calculate conductivity values and tensors from diffusion weighted imaging.
%
% Example inputs:
% dti_folder = '../data/CMRR_data_DBS_PSO/PD055_v2/dtifit';
% brain_seg_filename = '../data/CMRR_data_DBS_PSO/PD055_v2/PD055_FAST_seg_DWI_space.nii.gz';
%
% image_extent = [136,136,74];
% transform_to_DWI_space = reshape([1.5 0.000949223 0.000626029 0 0.00106894 -1.45873 -0.349425 0 0.000387685 0.349425 -1.45873 0 -98.3014 61.4612 103.401 1],4,4);
%
% saveFolder = 'PD055_STN_R/COMSOL_Models/ConductivityMaps';
%
% isSave = 1;
%
% version 2.1: optimized a little to run faster
%       * use load_nifti_v2 and send .gz files directly
%
% version 4: Skips version 3 (which incorporates partial skull) and builds
% directly off of version 2.1.  Allows conductivity to be input directly
% along with frequency. Otherwise, uses 6 kHz by default
%       * also fixes bug of relative permittivity for CSF (i.e. it should
%       be set to 109, or at least 1, but not 0 since it is supposed to be
%       a factor relative to free space)
%   
% version 5: adapts to slicer - focused modeling, removed Amira hacks 


function generate_conductivity_tensors_v4(dti_folder, brain_seg_filename, image_extent, transform_to_DWI_space, saveFolder, isSave, isotropic_electric_properties)

if nargin < 6
    warning('Warning: No isotropic electric properties specified.  Using those for 6 kHz by default...');
    freq = 6000;
    
    gray_sigma_iso = 0.11094;
    white_sigma_iso = 0.067393;
    
    gray_epsilon_iso = 35766;
    white_epsilon_iso = 18321;
else 
    freq = isotropic_electric_properties.freq;
    
    gray_sigma_iso = isotropic_electric_properties.gray_sigma_iso;
    white_sigma_iso = isotropic_electric_properties.white_sigma_iso;
    
    gray_epsilon_iso = isotropic_electric_properties.gray_epsilon_iso;
    white_epsilon_iso = isotropic_electric_properties.white_epsilon_iso;
end

csf_sigma_iso = 2;
csf_epsilon_iso = 109;


if ~exist(saveFolder,'dir')
    mkdir(saveFolder) % create folder to save files
end

% Load DTI eigenvalues and eigenvectors
L1 = load_nifti_v2([dti_folder filesep 'dti_L1.nii.gz']); % 1st eigenvalue
V1 = load_nifti_v2([dti_folder filesep 'dti_V1.nii.gz']); % 1st eigenvalue
L2 = load_nifti_v2([dti_folder filesep 'dti_L2.nii.gz']); % 2nd eigenvalue
V2 = load_nifti_v2([dti_folder filesep 'dti_V2.nii.gz']); % 2nd eigenvalue
L3 = load_nifti_v2([dti_folder filesep 'dti_L3.nii.gz']); % 3rd eigenvalue
V3 = load_nifti_v2([dti_folder filesep 'dti_V3.nii.gz']); % 3rd eigenvalue

nEl = numel(L1.vol);

[x_voxel_space,y_voxel_space,z_voxel_space] = ndgrid(0:V1.dim(2)-1 , 0:V1.dim(3)-1 , 0:V1.dim(4)-1);
x_voxel_space = x_voxel_space(:)';
y_voxel_space = y_voxel_space(:)';
z_voxel_space = z_voxel_space(:)';%(image_extent(3)-1)-z_voxel_space(:)'; % flip z axis since FSL and Amira read in different order; necessary for conversion

% Converts voxel indices into DWI space coordinates as viewed in Amira
x_DWI_space = ...
    transform_to_DWI_space(1,1)*x_voxel_space + transform_to_DWI_space(1,2)*y_voxel_space + transform_to_DWI_space(1,3)*z_voxel_space + transform_to_DWI_space(1,4);
y_DWI_space = ...
    transform_to_DWI_space(2,1)*x_voxel_space + transform_to_DWI_space(2,2)*y_voxel_space + transform_to_DWI_space(2,3)*z_voxel_space + transform_to_DWI_space(2,4);
z_DWI_space = ...
    transform_to_DWI_space(3,1)*x_voxel_space + transform_to_DWI_space(3,2)*y_voxel_space + transform_to_DWI_space(3,3)*z_voxel_space + transform_to_DWI_space(3,4);

xyz = [x_DWI_space;y_DWI_space;z_DWI_space];
clearvars x_DWI_space y_DWI_space z_DWI_space

L(1,1,:) = reshape(L1.vol,1,[]);
L(2,2,:) = reshape(L2.vol,1,[]);
L(3,3,:) = reshape(L3.vol,1,[]);
clearvars L1 L2 L3

V(1:3,1,:)=reshape(V1.vol,[],3)';
V(1:3,2,:)=reshape(V2.vol,[],3)';
V(1:3,3,:)=reshape(V3.vol,[],3)';
clearvars V1 V2 V3

% Flag voxels where all eigenvalues are zero or one is negative
invalid_z = L(1,1,:) == 0 & L(2,2,:) == 0 & L(3,3,:) == 0;
invalid_n = L(1,1,:) <= 0 | L(2,2,:) <= 0 | L(3,3,:) <= 0;
valid = squeeze(~(invalid_n | invalid_z));

tissue_type = load_nifti_v2(brain_seg_filename);
tissueType = reshape(tissue_type.vol,1,[]);

% gunzip(['../MR/' dtiG.Name '/FSL/T1_BrainOnly_seg_2diff-ventricles_mask.nii.gz']);
% ventricles_mask = load_nifti(['../MR/' dtiG.Name '/FSL/T1_BrainOnly_seg_2diff-ventricles_mask.nii']);
% ventriclesMask = reshape(ventricles_mask.vol,1,[]);

% Flag by tissue type and brain region in dark brain regions
% that can get missclassified as csf (SN, GP, and Cerebellar Nuclei)
csf_corrected = (tissueType == 1); % & ventriclesMask;
gray_corrected = (tissueType == 2);
white_corrected = tissueType == 3;
outsideBrain = tissueType == 0; %*** testing
valid = valid & ~outsideBrain';


% % % Tuch = 0.844e3; % Tuch et al (2001) scaling factor
% % % se_sigma = Tuch*L;

% % % %%
% % % figure
% % % h=histogram(se_sigma(1,1,valid));
% % %
% % % %%
% % % figure, hold on
% % % histogram(se_sigma(1,1,valid & gray_corrected'));
% % % histogram(se_sigma(1,1,valid & white_corrected'));
% % % histogram(se_sigma(1,1,valid & csf_corrected'));
% % %
% % % legend('gray','white','csf')

%%

% % % % Scaled eigenvalue (Tuch)
% % % Tuch = 0.844e3; % Tuch et al (2001) scaling factor
% % % se_sigma = Tuch*L;
% % % se_SIGMA = zeros(nEl,3,3);
% % % for ind = 1:nEl
% % %     se_SIGMA(ind,:,:) = V(:,:,ind) * se_sigma(:,:,ind) * transpose(V(:,:,ind));
% % % end

% Calculate conductivity values


epsilon = ones(nEl,1)*gray_epsilon_iso; % initalize all voxels as gray
epsilon(white_corrected) = white_epsilon_iso;
epsilon(csf_corrected) = csf_epsilon_iso;

% % % % Volume constraint (Wolters)
% % % vc_gray = eye(3)*gray_sigma_iso;
% % % vc_gray = repmat(vc_gray,[1 1 nEl]);
% % % 
% % % vc_csf = eye(3)*csf_sigma_iso;
% % % vc_csf = repmat(vc_csf,[1 1 nEl]);
% % % 
% % % 
% % % vc_white_perp = white_sigma_iso/nthroot(9,3);
% % % vc_white_para = vc_white_perp * 9;
% % % vc_white = [vc_white_para 0 0 ; 0 vc_white_perp 0 ; 0 0 vc_white_perp];
% % % vc_white = repmat(vc_white,[1 1 nEl]);
% % % 
% % % vc_sigma = vc_gray; % initalize all voxels as gray
% % % vc_sigma(:,:,csf_corrected) = vc_csf(:,:,csf_corrected);
% % % vc_sigma(:,:,white_corrected) = vc_white(:,:,white_corrected);


% Normalized volume constraint (Schmidt) sigma_i = d_i*(sigma_iso/nthroot(L1*L2*L3,3)
denom = squeeze(nthroot(L(1,1,:).*L(2,2,:).*L(3,3,:),3));
valid_nv = denom > 0; % flag voxels with crazy values

gray_norm_factor = gray_sigma_iso./denom;
white_norm_factor = white_sigma_iso./denom;

nv_csf = eye(3)*csf_sigma_iso;
nv_csf = repmat(nv_csf,[1 1 nEl]);

nv_gray = zeros(3,3,nEl);
nv_white = zeros(3,3,nEl);

% % % % Edits by Edgar Pe???a, 2/7/17 (part 1)
% % % % Isotropic tensor
% % % nv_gray_iso = zeros(3,3,nEl);
% % % nv_white_iso = zeros(3,3,nEl);
% % % % End Edits by Edgar Pe???a, 2/7/17 (part 1)

for iEl = 1:nEl
    
% % %     % Edits by Edgar Pe???a, 2/7/17 (part 2)
% % %     nv_gray_iso(:,:,iEl) = eye(3)*gray_sigma_iso;
% % %     nv_white_iso(:,:,iEl) = eye(3)*white_sigma_iso;
% % %     % End Edits by Edgar Pe???a, 2/7/17 (part 2)
    
    nv_gray(:,:,iEl) = L(:,:,iEl) * gray_norm_factor(iEl);
    nv_white(:,:,iEl) = L(:,:,iEl) * white_norm_factor(iEl);
end
clearvars L

% % % % Edits by Edgar Pe???a, 2/7/17 (part 2)
% % % nv_sigma_iso = nv_gray_iso;
% % % nv_sigma_iso(:,:,csf_corrected) = nv_csf(:,:,csf_corrected);
% % % nv_sigma_iso(:,:,white_corrected) = nv_white_iso(:,:,white_corrected);
% % % % End Edits by Edgar Pe???a, 2/7/17 (part 2)

nv_sigma = nv_gray; % initalize all voxels as gray
nv_sigma(:,:,csf_corrected) = nv_csf(:,:,csf_corrected);
nv_sigma(:,:,white_corrected) = nv_white(:,:,white_corrected);
clearvars nv_gray nv_csf nv_white

% calculate conductivity tensors
% vc_SIGMA = zeros(nEl,3,3);
nv_SIGMA = zeros(nEl,3,3);

% Edits by Edgar Pe???a, 2/7/17 (part 3)
% nv_SIGMA_iso = zeros(nEl,3,3);
% End Edits by Edgar Pe???a, 2/7/17 (part 3)
for iEl = 1:nEl
%     vc_SIGMA(iEl,:,:) = V(:,:,iEl) * vc_sigma(:,:,iEl) * transpose(V(:,:,iEl));
    nv_SIGMA(iEl,:,:) = V(:,:,iEl) * nv_sigma(:,:,iEl) * transpose(V(:,:,iEl));
    
    % Edits by Edgar Pe???a, 2/7/17 (part 4)
%     nv_SIGMA_iso(iEl,:,:) = V(:,:,iEl) * nv_sigma_iso(:,:,iEl) * transpose(V(:,:,iEl));
    % End Edits by Edgar Pe???a, 2/7/17 (part 4)
end
clearvars V nv_sigma

% vc = [xyz' vc_SIGMA(:,1,1) vc_SIGMA(:,1,2) vc_SIGMA(:,1,3) vc_SIGMA(:,2,2) vc_SIGMA(:,2,3) vc_SIGMA(:,3,3) epsilon];
nv = [xyz' nv_SIGMA(:,1,1) nv_SIGMA(:,1,2) nv_SIGMA(:,1,3) nv_SIGMA(:,2,2) nv_SIGMA(:,2,3) nv_SIGMA(:,3,3) epsilon];
% se = [xyz' se_SIGMA(:,1,1) se_SIGMA(:,1,2) se_SIGMA(:,1,3) se_SIGMA(:,2,2) se_SIGMA(:,2,3) se_SIGMA(:,3,3) epsilon];


% vcSorted = sortrows(vc(valid,:));
nvSorted = sortrows(nv(valid & valid_nv,:));
% seSorted = sortrows(se(valid,:));

if (isSave == 1)
    % save the resulting data file
%     save([saveFolder '/vc_' num2str(freq) '.txt'], 'vcSorted', '-ASCII');
    save([saveFolder '/nv_' num2str(round(freq)) '.txt'], 'nvSorted', '-ASCII');
%     save([saveFolder '/se_' num2str(freq) '.txt'], 'seSorted', '-ASCII');
else
    fprintf('Warning: isSave = 0\nResults will not be saved...\n');
end

