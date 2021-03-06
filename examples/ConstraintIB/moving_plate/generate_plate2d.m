%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2016 - 2019 by the IBAMR developers
%% All rights reserved.
%%
%% This file is part of IBAMR.
%%
%% IBAMR is free software and is distributed under the 3-clause BSD
%% license. The full text of the license can be found in the file
%% COPYRIGHT at the top level directory of IBAMR.
%%
%% ---------------------------------------------------------------------

clear all;
clc;

N = 1024;
H = 22;
dy = H/N;
b = 1.0;
Ny = ceil(b/dy);
Xcom = 0.0;

X = zeros(Ny,1);
Y = zeros(Ny,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Write out the vertex information
vertex_fid = fopen('plate2d.vertex', 'w');

% first line is the number of vertices in the file
fprintf(vertex_fid, '%d\n', Ny);

% remaining lines are the initial coordinates of each vertex
for k = 1:Ny

    X(k) = Xcom;
    Y(k) = -0.5 + (k-1)*dy;
    fprintf(vertex_fid, '%1.16e %1.16e\n', X(k), Y(k));
end %for
  

fclose(vertex_fid);

plot(X,Y,'o');

