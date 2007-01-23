%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters

NFINEST = 256;
num_nodes = 61*NFINEST/256;
X_center = [1.85 4.0];
radius = 0.15;

stiffness = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the vertex information
vertex_fid = fopen(["cylinder2d_" num2str(NFINEST) ".vertex"], "w");

% first line is number of vertices in the file
fprintf(vertex_fid, "%d\n", num_nodes);

% remaining lines are the initial coordinates of each vertex
for l = 0:num_nodes-1
  theta = 2.0*pi*l/num_nodes;
  X(1) = X_center(1) + radius*cos(theta);
  X(2) = X_center(2) + radius*sin(theta);
  fprintf(vertex_fid, "%1.16e %1.16e\n", X(1), X(2));
end %for

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 2: Write out the link information (including connectivity and
% material parameters).
edge_fid = fopen(["cylinder2d_" num2str(NFINEST) ".edge"], "w");

% first line is the number of edges
fprintf(edge_fid, "%d\n", num_nodes);

% remaining lines are the edges in the mesh
for l = 0:num_nodes-1
  current_idx = l;
  next_idx = current_idx+1;
  if (l == num_nodes-1)
    next_idx = 0;
  end %if

  kappa = stiffness*num_nodes;  % scale by 1/ds = num_nodes
  rest_length = 0.0;            % resting length of link

  fprintf(edge_fid, "%6d %6d %1.16e %1.16e\n", current_idx, next_idx, ...
          kappa, rest_length);
end %for

fclose(edge_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
