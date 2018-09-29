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