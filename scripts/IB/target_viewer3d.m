%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2008 - 2020 by the IBAMR developers
%% All rights reserved.
%%
%% This file is part of IBAMR.
%%
%% IBAMR is free software and is distributed under the 3-clause BSD
%% license. The full text of the license can be found in the file
%% COPYRIGHT at the top level directory of IBAMR.
%%
%% ---------------------------------------------------------------------

function target_viewer3d(structure_name)
  vertex_filename = [structure_name '.vertex'];
  target_filename = [structure_name '.target'];
  fprintf('reading vertex file named: %s\n', vertex_filename);
  fprintf('reading target point file named: %s\n', target_filename);

  fid = fopen(vertex_filename,'r');
  s = fgetl(fid);
  num_vertices = sscanf(s,'%d');
  X = zeros(num_vertices,1);
  Y = zeros(num_vertices,1);
  Z = zeros(num_vertices,1);
  for k = 1:num_vertices
    s = fgetl(fid);
    vals = sscanf(s,'%g %g %g');
    X(k) = vals(1);
    Y(k) = vals(2);
    Z(k) = vals(3);
  end

  fid = fopen(target_filename,'r');
  s = fgetl(fid);
  num_targets = sscanf(s,'%d');
  X_target = zeros(1,num_targets);
  for k = 1:num_targets
    s = fgetl(fid);
    val = sscanf(s,'%d');
    idx = val(1)+1;
    X_target(1,k) = X(idx);
    Y_target(1,k) = Y(idx);
    Z_target(1,k) = Z(idx);
  end

  hold on
  plot3(X_target,Y_target,Z_target,'k.')
  axis equal
  axis tight
  hold off