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

g = 9.81;
T = 5;
h = 10.0;

L_shallow = T*sqrt(g*h);

fun = @(L) L^2 - T^2 * g * (L/(2*pi))*tanh(2*pi*h/L);

L = fzero(fun, L_shallow)

%%
L = 2;
h = 0.5;
g = 9.81;
T = sqrt( L^2/(g * (L/(2*pi))*tanh(2*pi*h/L)))