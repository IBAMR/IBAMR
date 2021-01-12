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

PATH = './';

WallPath = strcat(PATH,'wedge3d.vertex');
%%
% Domain parameters
Lx = 12.0;   Ly = 2.0;     Lz = 3.0;
N = 50*2;
Nx = 6*N;   Ny = N; Nz = 3*N/2;

% Dimensional parameters
dx = Lx/Nx; dy = Ly/Ny; dz = Lz/Nz;
Length = 1.2; Width = 1.2;
Length_From_Left = Lx/2.0;
Distance_From_Edge = Ly/2.0;
Height_Above_Bottom = 2.3;
angle = 25; %degrees
x0 = Length_From_Left;
z0 = Height_Above_Bottom;

% From above, calculate the height
Height = Length/2*tan(25*pi/180);
m = Height/(Length/2); % slope of the xz line
br = z0-m*x0;
bl = z0+m*x0;

NumPtsX = ceil(Length/dx) + 1
NumPtsY = ceil(Width/dy) + 1
NumPtsZ = ceil(Height/dz) + 1

lag_pts = 0;

% Generate wedge shape
for k = 1:NumPtsZ
    z = Height_Above_Bottom + (k-1)*dz;
    for j = 1:NumPtsY
        y = Distance_From_Edge - Width/2 + (j-1)*dy;
        for i = 1:NumPtsX
            x = Length_From_Left-Length/2 + (i-1)*dx;
            Zr = m*x+br;
            Zl = -m*x+bl;
            if z >= Zr && z >= Zl
              lag_pts = lag_pts+1;
              LagX(lag_pts) = x;
              LagY(lag_pts) = y;
              LagZ(lag_pts) = z;
            end
        end
    end
end


figure;
plot3(LagX,LagY,LagZ,'.')
view(0,0)
axis equal;

%% write it in file

lag_pts = length(LagY);
fid = fopen(WallPath,'wt');
fprintf(fid,'%d\n',lag_pts);

for i = 1:length(LagX)
    fprintf(fid,'%12.7E\t\t%12.7E\t\t%12.7E\n',LagX(i), LagY(i),LagZ(i));
end

fclose(fid);
