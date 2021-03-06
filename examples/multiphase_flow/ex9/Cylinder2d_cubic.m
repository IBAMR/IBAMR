%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2019 - 2019 by the IBAMR developers
%% All rights reserved.
%%
%% This file is part of IBAMR.
%%
%% IBAMR is free software and is distributed under the 3-clause BSD
%% license. The full text of the license can be found in the file
%% COPYRIGHT at the top level directory of IBAMR.
%%
%% ---------------------------------------------------------------------

%% Filename : Cylinder2d_cubic.m

clf
clear all;
clc;
%% Mesh parameters.
radius = 0.5;
D = 2.0*radius
GAS_LS_INIT = 0.0
Lx = 5.0; Ly = 5.0;
Nx = 200; Ny = 200;
dx = Lx/Nx
dy = Ly/Ny


%% disk parameters.
h = dx;
X_shift = Lx/2; Y_shift = Ly/2;  
%h/16 to put marker away from grid lines on second level
%The markers are symmetric on single level (by taking odd #cells)
X_com = 0.0;
Y_com = 0.0;

%% write out the points.

counter = 0;
Npts = floor(radius/h);

%1st quadrant
idx = 0; 
for j = 1:Npts
    y = Y_com + (j)*h;
    for i = 1:Npts
        x = X_com + (i)*h;
        
        dist2 = (x-X_com)^2 + (y-Y_com)^2;
        if (dist2 <= radius^2)
            idx = idx+1;
            X_first(idx) = x; 
            Y_first(idx) = y; 
        end
    end
end
quad_pts = idx;

%2nd quadrant
idx = 0;
for k = 1:quad_pts
            x = X_first(k);
            y = Y_first(k);
          
            idx = idx+1;
            X_second(idx) = -x; 
            Y_second(idx) = y; 
end

%3rd quadrant
idx  = 0;
for k = 1:quad_pts
            x = X_first(k);
            y = Y_first(k);
          
            idx = idx+1;
            X_third(idx) = -x; 
            Y_third(idx) = -y; 
end

%4th quadrant
idx  = 0;
for k = 1:quad_pts
            x = X_first(k);
            y = Y_first(k);
          
            idx = idx+1;
            X_fourth(idx) = x; 
            Y_fourth(idx) = -y; 
end

%X-axis
idx = 0;
for k = -Npts+1:Npts-1
            x = X_com + k*h;
            y = Y_com;
          
            idx = idx+1;
            X_xaxis(idx) = x; 
            Y_xaxis(idx) = y; 
end

%Y-axis
idx = 0;
for k = -Npts+1:Npts-1
            x = X_com;
            y = Y_com + k*h;
          
            idx = idx+1;
            X_yaxis(idx) = x; 
            Y_yaxis(idx) = y; 
end

%%

X_array = [X_first,X_second,X_third,X_fourth,X_xaxis,X_yaxis];
Y_array = [Y_first,Y_second,Y_third,Y_fourth,Y_xaxis,Y_yaxis];

X_array = X_array + X_shift;
Y_array = Y_array + Y_shift;
    
%% plot the coordinates
clf
plot(X_first, Y_first, 'r.');
hold on;
plot(X_second, Y_second, 'k.');
hold on;
plot(X_third, Y_third, 'g.');
hold on;
plot(X_fourth, Y_fourth, 'c.');
hold on;
plot(X_xaxis, Y_xaxis, 'ro');
hold on;
plot(X_yaxis, Y_yaxis, 'ko');
hold on;
axis('equal')
hold off;


%% write the coordinates in the file
fid = fopen('cylinder2d.vertex','wt');
fprintf(fid,'%d\n', length(X_array));
for i = 1:length(X_array)
    fprintf(fid,'%f\t%f\n',X_array(i),Y_array(i));
end
fclose(fid); 

%% plot the disk

clf
plot(X_array, Y_array, 'r.');