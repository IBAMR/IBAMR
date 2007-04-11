function generate_mesh2d(NFINEST,generate_coarse_mesh)

% NOTE: NFINEST = 8 corresponds to a uniform grid spacing of h=1/128

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters

num_nodes = 50*NFINEST;
num_layers = 4*NFINEST;

width = 8.0/64.0;
taper_width = 3.0/256.0;

alpha = 0.2;
beta  = 0.275;

stiffness = 1.0/num_layers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Write out the vertex information
if (generate_coarse_mesh)
  vertex_fid = fopen(['shell2d_crse_' num2str(16*NFINEST) '.vertex'], 'w');
else
  vertex_fid = fopen(['shell2d_fine_' num2str(16*NFINEST) '.vertex'], 'w');
end %if

% determine the number of layers in the mesh
num_coarse_layers = 0;
num_fine_layers = 0;
for r = 0:num_layers
  delta_r = width*(r/num_layers - 0.5);

  if (delta_r > -0.5*width+taper_width & ...
      delta_r < +0.5*width-taper_width)
    num_coarse_layers = num_coarse_layers+1;
  end %if

  if (delta_r < -0.5*width+2.0*taper_width | ...
      delta_r > +0.5*width-2.0*taper_width)
    num_fine_layers = num_fine_layers+1;
  end %if

end %for

if (generate_coarse_mesh)
  num_mesh_layers = num_coarse_layers;
else
  num_mesh_layers = num_fine_layers;
end %if

% first line is the number of vertices in the file
fprintf(vertex_fid, '%d\n', num_mesh_layers*num_nodes);

% remaining lines are the initial coordinates of each vertex
for r = 0:num_layers
  delta_r = width*(r/num_layers - 0.5);

  use_layer = false;
  if (generate_coarse_mesh & ...
      (delta_r > -0.5*width+taper_width & ...
       delta_r < +0.5*width-taper_width))
    use_layer = true;
  elseif (~generate_coarse_mesh &
          (delta_r < -0.5*width+2.0*taper_width | ...
           delta_r > +0.5*width-2.0*taper_width))
    use_layer = true;
  end %if

  if (use_layer)
    for l = 0:num_nodes-1
      theta = 2.0*pi*l/num_nodes;
      X(1) = 0.5 + (alpha+delta_r)*cos(theta);
      X(2) = 0.5 + (beta +delta_r)*sin(theta);
      fprintf(vertex_fid, '%1.16e %1.16e\n', X(1), X(2));
    end %for
  end %if
end %for

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 2: Write out the link information (including connectivity and
% material parameters).
if (generate_coarse_mesh)
  spring_fid = fopen(['shell2d_crse_' num2str(16*NFINEST) '.spring'], 'w');
else
  spring_fid = fopen(['shell2d_fine_' num2str(16*NFINEST) '.spring'], 'w');
end %if

% first line is the number of edges in the file
fprintf(spring_fid, '%d\n', num_mesh_layers*num_nodes);

offset = 0;

% remaining lines are the edges in the mesh
for r = 0:num_layers
  delta_r = width*(r/num_layers - 0.5);

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

  if (~generate_coarse_mesh)
    r_fac = 1.0-r_fac;
  end %if

  use_layer = false;
  if (generate_coarse_mesh & ...
      (delta_r > -0.5*width+taper_width & ...
       delta_r < +0.5*width-taper_width))
    use_layer = true;
  elseif (~generate_coarse_mesh &
          (delta_r < -0.5*width+2.0*taper_width | ...
           delta_r > +0.5*width-2.0*taper_width))
    use_layer = true;
  end %if

  if (use_layer)
    for l = 0:num_nodes-1
      current_idx = l + offset;
      next_idx = current_idx+1;
      if (l == num_nodes-1)
        next_idx = 0 + offset;
      end %if

      kappa = r_fac*stiffness*num_nodes;  % scale by 1/ds = num_nodes
      rest_length = 0.0;                  % resting length of link

      fprintf(spring_fid, '%6d %6d %1.16e %1.16e\n', current_idx, next_idx, ...
              kappa, rest_length);
    end %for

    offset = offset + num_nodes;
  end %if
end %for

fclose(spring_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
