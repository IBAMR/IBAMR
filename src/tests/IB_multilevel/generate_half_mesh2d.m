%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters

NFINEST = 4;  % NFINEST = 4 corresponds to a uniform grid spacing of h=1/64
num_nodes = 25*NFINEST+1;
num_layers = 2*NFINEST;

width = 4.0/64.0;
alpha = 0.2;
beta  = 0.25;

tapered_stiffness = false;
stiffness = 1.0/num_layers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Write out the vertex information
vertex_fid = fopen("shell2d_right_64.vertices", "w");

% first line is number of vertices in the file
fprintf(vertex_fid, "%d\n", num_layers*num_nodes);

% remaining lines are the initial coordinates of each vertex
for r = 0:num_layers-1
  for l = 0:num_nodes-1
    theta = pi*l/(num_nodes-1) - pi/2;
    delta_r = width*((r+0.5)/num_layers - 0.5);
    X(1) = 0.5 + (alpha+delta_r)*cos(theta);
    X(2) = 0.5 + (beta +delta_r)*sin(theta);
    fprintf(vertex_fid, "%1.16e %1.16e\n", X(1), X(2));
  end %for
end %for

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 2: Write out the link information (including connectivity and
% material parameters).
edge_fid = fopen("shell2d_right_64.edges", "w");

% first line is the base index; here we use 0-based indexing
fprintf(edge_fid, "%d\n", 0);

% second line is the number of edges
fprintf(edge_fid, "%d\n", num_layers*(num_nodes-1));

% remaining lines are the edges in the mesh
for r = 0:num_layers-1
  if (tapered_stiffness)
    theta = -0.5 + (r+0.5)/num_layers;
    r_fac = 1.0 + cos(2.0*pi*theta);
  else
    r_fac = 1.0;
  end %if

  for l = 0:num_nodes-2
    current_idx = l + r*num_nodes;
    next_idx = current_idx+1;

    kappa = r_fac*stiffness*2*(num_nodes-1);  % scale by 1/ds = num_nodes
    rest_length = 0.0;                        % resting length of link

    fprintf(edge_fid, "%6d %6d %1.16e %1.16e\n", current_idx, next_idx, ...
            kappa, rest_length);
  end %for
end %for

fclose(edge_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
