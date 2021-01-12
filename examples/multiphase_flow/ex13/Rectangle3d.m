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
  
  WallPath = strcat(PATH,'rectangle3d.vertex');
  %%
  % Domain parameters
  Lx = 3.22; Ly = 1.0; Lz = 1.0;
  Nx = 161;  Ny = 50; Nz = 50;
  
  % Dimensional parameters
  dx = Lx/Nx; dy = Ly/Ny; dz = Lz/Nz;
  Length = 0.16;
  Width  = 0.4;
  Height = 0.16;
  Xcom   = 0.74;
  Ycom   = 0.3 + Width/2;
  Zcom   = Height/2.0;
  
  NumPtsX = ceil(Length/dx) + 1
  NumPtsY = ceil(Width/dy) + 1
  NumPtsZ = ceil(Height/dz) + 1
  
  lag_pts = 0;
  
  % Generate Rectangular Prism
  for k = 1:NumPtsZ
    z = Zcom - Height/2 + (k-1)*dz;
      for j = 1:NumPtsY
          y = Ycom - Width/2 + (j-1)*dy;
          for i = 1:NumPtsX
                x = Xcom-Length/2 + (i-1)*dy;
                lag_pts = lag_pts+1;
                LagX(lag_pts) = x;
                LagY(lag_pts) = y;
                LagZ(lag_pts) = z;
          end
      end
  end
  

  figure;
  plot3(LagX,LagY, LagZ,'.')
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
