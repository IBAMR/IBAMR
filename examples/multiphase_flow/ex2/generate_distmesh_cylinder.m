clear all;
clc;

% Add path to distmesh directory
%http://persson.berkeley.edu/distmesh/
%addpath('/path/to/distmesh');

R           = 0.055;
D           = 2.0*R;
LENGTH      = 20.0*D;
Nx          = 880;
H           = LENGTH/Nx;

x0 = [0, 0];  
r = R;
h  = H;
xlower = -r + x0(1);
ylower = -r + x0(2);
xupper = r + x0(1);
yupper = r + x0(2);

% Define implicit surface.
fd = @(p) sqrt( sum( (p - x0).^2,2 ) ) - r;
[p,t] = distmesh2d(fd, @huniform, h , [xlower, ylower; xupper, yupper], []);

% Format the arrays that libMesh expects.
domain_id = 0;
t = [t,domain_id *ones(length(t),1)];
p = p';
t = t';

filename = 'cylinder.xda';
fid = fopen(filename, 'w'); 
fprintf(fid, '%d %d \n', length(p), length(t)); 
fprintf(fid, '%f %f \n', p); 
fprintf(fid, '%d %d %d %d \n', t); 
fclose(fid);




