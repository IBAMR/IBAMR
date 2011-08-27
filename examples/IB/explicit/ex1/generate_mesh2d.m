%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters

N = 64;
beta  = 0.25;
alpha = 0.25^2/beta;

A = pi*alpha*beta;      % area of ellipse
R = sqrt(A/pi);         % radius of disc with equivalent area as the ellipse
perim = 2*pi*R;         % perimeter of the equivalent disc
dx = 1/N;
num_nodes = ceil(perim/(dx/3)/4)*4;
ds = 1.0/num_nodes;

kappa = 1.0/perim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Write out the vertex information
vertex_fid = fopen(['curve2d_' num2str(N) '.vertex'], 'w');

% first line is the number of vertices in the file
fprintf(vertex_fid, '%d\n', num_nodes);

% remaining lines are the initial coordinates of each vertex
for l = 0:num_nodes-1
  theta = 2.0*pi*l/num_nodes;
  X(1) = 0.5 + (alpha)*cos(theta);
  X(2) = 0.5 + (beta )*sin(theta);
  fprintf(vertex_fid, '%1.16e %1.16e\n', X(1), X(2));
end %for

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 2: Write out the link information (including connectivity and
% material parameters).
spring_fid = fopen(['curve2d_' num2str(N) '.spring'], 'w');

% first line is the number of edges in the file
fprintf(spring_fid, '%d\n', num_nodes);

% remaining lines are the edges in the mesh
for l = 0:num_nodes-1
  current_idx = l;
  next_idx    = current_idx+1;
  if (l == num_nodes-1)
    next_idx = 0;
  end %if
  
  K = kappa/ds;
  rest_length = 0.0;
  
  fprintf(spring_fid, '%6d %6d %1.16e %1.16e\n', current_idx, next_idx, ...
          K, rest_length);
end %for

fclose(spring_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
