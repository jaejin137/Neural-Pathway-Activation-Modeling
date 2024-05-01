% COMSOL_to_DWI_transform = ...
%     [-.43 -.84 .33 9.35;
%     .87 -.28 .41 23.86;
%     -.25 .46 .85 5.77;
%     0 0 0 1];
% COMSOL_to_DWI_transform = ...
%    [1 0 0 0; 0 0.98 0 0; 0 0 1 0; 0 0 0 1]
    COMSOL_to_DWI_transform = ...
    [0.673416 0.723739 0.0393856 -2.19027;
    -0.680111 0.586049 0.424406 7.56046;
    0.289759 -0.306459 0.904615 -0.196531 ]

format compact

rot_matrix = COMSOL_to_DWI_transform (1:3,1:3)
rot_vector = vrrotmat2vec(rot_matrix)
%rot_vector = rotmat2vec3d(rot_matrix)
rx = rot_vector(1)
ry = rot_vector(2)
rz = rot_vector(3)
angle = rad2deg(rot_vector(4))
sx = COMSOL_to_DWI_transform(1,4) 
sy = COMSOL_to_DWI_transform(2,4) 
sz = COMSOL_to_DWI_transform(3,4) 