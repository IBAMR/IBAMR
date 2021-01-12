%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2019 - 2020 by the IBAMR developers
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
%% Mesh parameters.
depth = 0.4; 
Lx = 68.0 * depth; Ly = 1.3 * depth;
Nx = 1768; Ny = 264;
dx = Lx/Nx; dy = Ly/Ny;

%% Trapzoid parameters.
Lt = 27.5 * depth;
ht = 0.75 * depth;
num_pts_x = ceil(Lt/dx);
num_pts_y = ceil(ht/dy);

%% Bottom left and bottom right corner
x1 = 15.0 * depth;
y1 = 0.0;
x2 = x1 + Lt;
y2 = 0.0;

%% Compute the parameters for the "cut out" lines
l1 = 15 * depth;
m1 = ht / l1;
b1 = y1 - m1 * x1;

l2 = 7.5 * depth;
m2 = -ht / l2;
b2 = y2 - m2 * x2;

%% write out the points.
idx = 0;
for i = 1:num_pts_x
    x = x1 + ( (i-1)*dx );
    
    for j = 1:num_pts_y
        y = y1 + ( (j-1)*dy );
        
        if( y <= m1 * x + b1 && y <= m2 * x + b2)
            idx = idx+1;
            X_array(idx) = x; Y_array(idx) = y;
        end
    end
end

XCOM = sum(X_array)/length(X_array);
YCOM = sum(Y_array)/length(Y_array);

X_array = X_array;
Y_array = Y_array;
%% plot the trap
plot(X_array, Y_array, '.');
%axis 'equal'

%% write the coordinates in the file
fid = fopen('trapezoid2d.vertex','wt');
fprintf(fid,'%d\n', length(X_array));

for i = 1:length(X_array)
    fprintf(fid,'%f\t%f\n',X_array(i),Y_array(i));
end