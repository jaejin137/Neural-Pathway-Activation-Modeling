function in = is_inside(volume, q)

% Define the vertices and faces of a tetrahedron
% Define the vertices and faces of a cube
%v = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
%f = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

% Create a structure with vertices and faces
%fv.vertices = v;
%fv.faces = f;

% Define the query point
%q = [0.5 0.5 0.5];

% Test if the point is inside the tetrahedron
in = inpolyhedron(fv, q);

% Display the result
% if in
%     disp('The point is inside the tetrahedron.')
% else
%     disp('The point is outside the tetrahedron.')
% end