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
  
  WallPath = strcat(PATH,'barge2d.vertex');
  %%
  % Domain parameters
  Lx = 5.0; Ly = 2.5;
  Nx = 500*2*2;  Ny = 250*2*2;
  
  % Dimensional parameters
  dx = Lx/Nx; dy = Ly/Ny;
  Length = 0.3;
  Width  = Length/3.0;
  Xcom   = 0.0;
  Ycom   = 0.0;
  theta  = 15; %degrees
  t      = theta*pi/180; % radians
  
  NumPtsX = ceil(Length/dx) + 1
  NumPtsY = ceil(Width/dy) + 1
  
  lag_pts = 0;
  
  % Generate Rectangular Prism
  for j = 1:NumPtsY
      y = Ycom - Width/2 + (j-1)*dy;
          for i = 1:NumPtsX
              x = Xcom-Length/2 + (i-1)*dy;
                lag_pts = lag_pts+1;
                LagX(lag_pts) = x;
                LagY(lag_pts) = y;
          end
  end
  
  % Rotate by theta
  Rx = LagX * cos(t) - LagY * sin(t);
  Ry = LagX * sin(t) + LagY * cos(t);
  LagX = Rx;
  LagY = Ry;


  figure;
  plot(LagX,LagY,'.')
  axis equal;
  
  %% write it in file
  
  lag_pts = length(LagY);
  fid = fopen(WallPath,'wt');
  fprintf(fid,'%d\n',lag_pts);
  
  for i = 1:length(LagX)
      fprintf(fid,'%12.7E\t\t%12.7E\n',LagX(i),LagY(i));
  end
  
  fclose(fid);
