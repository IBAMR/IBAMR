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

function vertex_viewer3d(structure_name)
  vertex_filename = [structure_name '.vertex'];
  fprintf('reading vertex file named: %s\n', vertex_filename);
  fid = fopen(vertex_filename,'r');
  s = fgetl(fid);
  num_vertices = sscanf(s,'%d');

  X = zeros(1,num_vertices);
  Y = zeros(1,num_vertices);
  Z = zeros(1,num_vertices);
  for k = 1:num_vertices
    s = fgetl(fid);
    vals = sscanf(s,'%g %g %g');
    X(1,k) = vals(1);
    Y(1,k) = vals(2);
    Z(1,k) = vals(3);
  end

  hold on
  plot3(X,Y,Z,'b.')
  axis equal
  axis tight
  hold off