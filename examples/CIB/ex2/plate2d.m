%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2015 - 2019 by the IBAMR developers
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

Lx = 32;
Ly = 22;
b = 1.0;
Nx = 256;
Ny = 256;

Dx = Lx / Nx;
Dy = Ly / Ny;
total_points = ceil(b/(1.25*Dy));

X_com = 10.0;
X = X_com*ones(total_points, 1);
Y = zeros(total_points, 1);

for j = 1: total_points
    Y(j,1) =  10.0 + (j-1)*Dy;
end

plot(X(:,1), Y(:,1), 'o');

%% Print out vertex files
vertex_fid = fopen('plate2d.vertex', 'w');
fprintf(vertex_fid, '%d\n', total_points);

for node = 1: total_points
     fprintf(vertex_fid, '%1.16e %1.16e\n', X(node,1), Y(node,1));
end 
