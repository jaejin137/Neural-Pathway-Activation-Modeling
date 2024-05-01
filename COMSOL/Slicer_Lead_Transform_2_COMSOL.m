COMSOL_to_DWI_transform = ...
    [-.43 -.84 .33 9.35;
    .87 -.28 .41 23.86;
    -.25 .46 .85 5.77;
    0 0 0 1];

format compact

rot_matrix = COMSOL_to_DWI_transform (1:3,1:3)
rot_vector = vrrotmat2vec(rot_matrix)
rx = rot_vector(1)
ry = rot_vector(2)
rz = rot_vector(3)
angle = rad2deg(rot_vector(4))
sx = COMSOL_to_DWI_transform(1,4) 
sy = COMSOL_to_DWI_transform(2,4) 
sz = COMSOL_to_DWI_transform(3,4) 