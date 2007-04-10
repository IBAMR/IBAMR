%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters
R = 0.5;
D = 2.0*R;
L_x = 20.0*D;
L_z = pi*D;
NFINEST_x = 192;
NFINEST_z = 32;
dx = L_x/NFINEST_x;
dz = L_z/NFINEST_z;

% choose num_nodes_r so that dr = L_z/num_nodes = dz (approximately)
num_nodes_r = NFINEST_z;
dr = dz;

% choose num_nodes_s so that ds = 2*PI*R/num_nodes = dx (approximately)
num_nodes_s = 0.1875*NFINEST_x;
ds = 2*pi*R/num_nodes_s;

X_center = [0.0 0.0 0.0];
stiffness = 1.0;

sprintf('dx = %f\n',dx);
sprintf('dr = %f\n',dr);
sprintf('ds = %f\n',ds);
sprintf('dz/dr = %f\n',dz/dr);
sprintf('dx/ds = %f\n',dx/ds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the vertex information
vertex_fid = fopen(['cylinder3d_' num2str(NFINEST_x) '.vertex'], 'w');

% first line is number of vertices in the file
fprintf(vertex_fid, '%d\n', num_nodes);

% remaining lines are the initial coordinates of each vertex
for k = 0:num_nodes_r-1
  for l = 0:num_nodes_s-1
    theta = 2.0*pi*l/num_nodes;
    X(1) = X_center(1) + R*cos(theta);
    X(2) = X_center(2) + R*sin(theta);
    X(3) = X_center(3) + (k+0.5)*dr - 0.5*L_z;
    fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', X(1), X(2), X(3));
  end %for
end %for

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0

% Step 2: Write out the link information (including connectivity and
% material parameters).
spring_fid = fopen(['cylinder3d_' num2str(NFINEST_x) '.spring'], 'w');

% first line is the number of edges
fprintf(spring_fid, '%d\n', num_nodes);

% remaining lines are the edges in the mesh
for l = 0:num_nodes-1
  current_idx = l;
  next_idx = current_idx+1;
  if (l == num_nodes-1)
    next_idx = 0;
  end %if

  kappa = stiffness/ds;  % scale by 1/ds
  rest_length = ds;      % resting length of link

  fprintf(spring_fid, '%6d %6d %1.16e %1.16e\n', current_idx, next_idx, ...
          kappa, rest_length);
end %for

fclose(spring_fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
