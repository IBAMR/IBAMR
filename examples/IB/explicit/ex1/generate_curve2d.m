%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2014 - 2019 by the IBAMR developers
%% All rights reserved.
%%
%% This file is part of IBAMR.
%%
%% IBAMR is free software and is distributed under the 3-clause BSD
%% license. The full text of the license can be found in the file
%% COPYRIGHT at the top level directory of IBAMR.
%%
%% ---------------------------------------------------------------------

function generate_curve2d(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters
beta  = 0.35;
alpha = 0.25^2/beta;

A = pi*alpha*beta;      % area of ellipse
R = sqrt(A/pi);         % radius of disc with equivalent area as the ellipse
perim = 2*pi*R;         % perimeter of the equivalent disc
dx = 1/N;
dx_64 = 1/64;
num_node_circum = (dx_64/dx)*ceil(perim/(dx_64/3)/4)*4;
ds = 2.0*pi*R/num_node_circum;
num_node = num_node_circum;

kappa = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Write out the vertex information
vertex_fid = fopen(['curve2d_' num2str(N) '.vertex'], 'w');

% first line is the number of vertices in the file
fprintf(vertex_fid, '%d\n', num_node);

% remaining lines are the initial coordinates of each vertex
for l = 0:num_node_circum-1
  theta = 2.0*pi*l/num_node_circum;
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
fprintf(spring_fid, '%d\n', num_node_circum);

% remaining lines are the edges in the mesh
for l = 0:num_node_circum-1
  current_idx = l;
  next_idx    = current_idx+1;
  if (l == num_node_circum-1)
    next_idx = 0;
  end %if
  
  K = kappa/ds;
  rest_length = 0.0;
  
  fprintf(spring_fid, '%6d %6d %1.16e %1.16e\n', current_idx, next_idx, ...
          K, rest_length);
end %for

fclose(spring_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
