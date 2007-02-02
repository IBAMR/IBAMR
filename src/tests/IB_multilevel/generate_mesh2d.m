%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters

NFINEST = 8;                                                               % NFINEST = 8 corresponds to a uniform grid spacing of h=1/128
num_nodes = 50*NFINEST;
num_layers = 2*NFINEST;

width = 4.0/64.0;
taper_width = 1.0/128.0;

alpha = 0.2;
beta  = 0.3;

coarse_level_stiffness = true;
stiffness = 1.0/num_layers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Write out the vertex information
if (coarse_level_stiffness)
  vertex_fid = fopen(["shell2d_crse_" num2str(16*NFINEST) ".vertex"], "w");
else
  vertex_fid = fopen(["shell2d_fine_" num2str(16*NFINEST) ".vertex"], "w");
end %if

% first line is the number of vertices in the file
fprintf(vertex_fid, "%d\n", num_layers*num_nodes);

% remaining lines are the initial coordinates of each vertex
for r = 0:num_layers-1
  delta_r = width*((r+0.5)/num_layers - 0.5);
  for l = 0:num_nodes-1
    theta = 2.0*pi*l/num_nodes;
    X(1) = 0.5 + (alpha+delta_r)*cos(theta);
    X(2) = 0.5 + (beta +delta_r)*sin(theta);
    fprintf(vertex_fid, "%1.16e %1.16e\n", X(1), X(2));
  end %for
end %for

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 2: Write out the link information (including connectivity and
% material parameters).
if (coarse_level_stiffness)
  edge_fid = fopen(["shell2d_crse_" num2str(16*NFINEST) ".edge"], "w");
else
  edge_fid = fopen(["shell2d_fine_" num2str(16*NFINEST) ".edge"], "w");
end %if

% first line is the number of edges in the file
fprintf(edge_fid, "%d\n", num_layers*num_nodes);

% remaining lines are the edges in the mesh
for r = 0:num_layers-1
  delta_r = width*((r+0.5)/num_layers - 0.5);

  if (delta_r < -0.5*width+taper_width)
    r_fac = 0.0;
  elseif (delta_r < -0.5*width+2.0*taper_width)
    r_fac = (delta_r+0.5*width-taper_width)/taper_width;
  elseif (delta_r < 0.5*width-2.0*taper_width)
    r_fac = 1.0;
  elseif (delta_r < 0.5*width-taper_width)
    r_fac = (0.5*width-taper_width-delta_r)/taper_width;
  else
    r_fac = 0.0;
  end %if

  if (~coarse_level_stiffness)
    r_fac = 1.0-r_fac;
  end %if

  for l = 0:num_nodes-1
    current_idx = l + r*num_nodes;
    next_idx = current_idx+1;
    if (l == num_nodes-1)
      next_idx = 0 + r*num_nodes;
    end %if

    kappa = r_fac*stiffness*num_nodes;  % scale by 1/ds = num_nodes
    rest_length = 0.0;                  % resting length of link

    fprintf(edge_fid, "%6d %6d %1.16e %1.16e\n", current_idx, next_idx, ...
            kappa, rest_length);
  end %for
end %for

fclose(edge_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
