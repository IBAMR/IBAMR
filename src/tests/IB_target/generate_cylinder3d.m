%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters
R = 0.5;
L = 16.0;
L_z = 4.0;
NFINEST = 128;
NFINEST_z = NFINEST/4;
dx = L/NFINEST;
dz = L_z/NFINEST_z;

% choose num_nodes_s so that ds = 2*PI*R/num_nodes = 0.5*dx (approximately)
num_nodes_s = 0.375*NFINEST;
ds = 2*pi*R/num_nodes_s;

% choose num_nodes_r so that dr = L_z/num_nodes_z = 0.5*dz (approximately)
num_nodes_r = 2.0*NFINEST_z;
dr = L_z/num_nodes_r;

X_center = [0.0 0.0 0.0];
stiffness = 1.0;

sprintf('dx = %f\n',dx);
sprintf('dr = %f\n',dr);
sprintf('ds = %f\n',ds);
sprintf('dz/dr = %f\n',dz/dr);
sprintf('dx/ds = %f\n',dx/ds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the vertex information
vertex_fid = fopen(['cylinder3d_' num2str(NFINEST) '.vertex'], 'w');

% first line is number of vertices in the file
fprintf(vertex_fid, '%d\n', num_nodes_r*num_nodes_s);

% remaining lines are the initial coordinates of each vertex
for k = 0:num_nodes_r-1
  for l = 0:num_nodes_s-1
    theta = 2.0*pi*l/num_nodes_s;
    X(1) = X_center(1) + R*cos(theta);
    X(2) = X_center(2) + R*sin(theta);
    X(3) = X_center(3) + (k+0.5)*dr - 0.5*L_z;
    fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', X(1), X(2), X(3));
  end %for
end %for

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
