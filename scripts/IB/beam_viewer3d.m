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

function beam_viewer3d(structure_name)
  vertex_filename = [structure_name '.vertex'];
  beam_filename = [structure_name '.beam'];
  fprintf('reading vertex file named: %s\n', vertex_filename);
  fprintf('reading beam file named: %s\n', beam_filename);

  fid = fopen(vertex_filename,'r');
  s = fgetl(fid);
  num_vertices = sscanf(s,'%d');
  X = zeros(num_vertices,1);
  Y = zeros(num_vertices,1);
  Z = zeros(num_vertices,1);
  for k = 1:num_vertices
    s = fgetl(fid);
    vals = sscanf(s,'%g');
    X(k) = vals(1);
    Y(k) = vals(2);
    Z(k) = vals(3);
  end

  fid = fopen(beam_filename,'r');
  s = fgetl(fid);
  num_beams = sscanf(s,'%d');
  X_beam = zeros(3,num_beams);
  Y_beam = zeros(3,num_beams);
  Z_beam = zeros(3,num_beams);
  for k = 1:num_beams
    s = fgetl(fid);
    vals = sscanf(s,'%d %d %d');
    idx1 = vals(1)+1;
    idx2 = vals(2)+1;
    idx3 = vals(3)+1;
    X_beam(1,k) = X(idx1);
    X_beam(2,k) = X(idx2);
    X_beam(3,k) = X(idx3);
    Y_beam(1,k) = Y(idx1);
    Y_beam(2,k) = Y(idx2);
    Y_beam(3,k) = Y(idx3);
    Z_beam(1,k) = Z(idx1);
    Z_beam(2,k) = Z(idx2);
    Z_beam(3,k) = Z(idx3);
  end

  hold on
  plot3(X_beam,Y_beam,Z_beam,'b')
  axis equal
  axis tight
  hold off