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

function spring_viewer3d(structure_name)
  vertex_filename = [structure_name '.vertex'];
  spring_filename = [structure_name '.spring'];
  fprintf('reading vertex file named: %s\n', vertex_filename);
  fprintf('reading spring file named: %s\n', spring_filename);

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

  fid = fopen(spring_filename,'r');
  s = fgetl(fid);
  num_springs = sscanf(s,'%d');
  X_spring = zeros(2,num_springs);
  Y_spring = zeros(2,num_springs);
  Z_spring = zeros(2,num_springs);
  for k = 1:num_springs
    s = fgetl(fid);
    vals = sscanf(s,'%d %d');
    idx1 = vals(1)+1;
    idx2 = vals(2)+1;
    X_spring(1,k) = X(idx1);
    X_spring(2,k) = X(idx2);
    Y_spring(1,k) = Y(idx1);
    Y_spring(2,k) = Y(idx2);
    Z_spring(1,k) = Z(idx1);
    Z_spring(2,k) = Z(idx2);
  end

  hold on
  plot3(X_spring,Y_spring,Z_spring,'b')
  axis equal
  axis tight
  hold off