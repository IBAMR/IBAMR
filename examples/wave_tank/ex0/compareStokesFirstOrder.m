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

clear all;

clc;

DEPTH       =  0.5;

AMPLITUDE   =  0.025/2.0;

TIME_PERIOD =  1.5;

OMEGA       =  2*pi/TIME_PERIOD;

WAVELENGTH  =  2.8264;               % From Stokes theory

WAVENUMBER  =  2.223;

TIME = linspace(0,60,1000);

%x    = 2.0;
x     = 0.065;

D = importdata('probe_0');
A = D.data


%phase = 0.81; %0.38;
phase = 0.0;

y = AMPLITUDE*cos(WAVENUMBER*x - OMEGA*(TIME - phase));



% If comparing from a cell above the still water depth, need to plot -phi + phi(0)
% where h is the uniform grid spacing
t = A(:,1);
eta = -A(:,2) + A(1,2);
plot(TIME,y,'k-', t, eta,'r-');

% Legend
alpha = 2.0;
xlabel('Time');
ylabel('Elevation');
legend('Analytical Stokes 1st', 'Numerical');
title(['x = ' num2str(x) ', alpha = ' num2str(alpha)]);