clear all;

clc;


dx = 0.01;

DEPTH       =  0.5;

AMPLITUDE   =  0.025/2.0;

TIME_PERIOD =  1.5;

OMEGA       =  2*pi/TIME_PERIOD;

WAVELENGTH  =  2.8264;               % From Stokes theory

WAVENUMBER  =  2.223;

TIME = linspace(0,60,1000);

%x    = 2.0;
x     = 3.165;

D = importdata('probe_1');
A = D.data;


%phase = 0.81; %0.38;
phase = 0.0;

H = 2*AMPLITUDE;
y = H/2*cos(WAVENUMBER*x - OMEGA*(TIME - phase)) + ... 
    pi*H^2/(8*WAVELENGTH) * cosh(WAVENUMBER*DEPTH) * ...
    (2+cosh(2*WAVENUMBER*DEPTH))/(sinh(WAVENUMBER*DEPTH)^3) * ... 
    cos(2*(WAVENUMBER*x - OMEGA*(TIME - phase)));



% If comparing from a cell above the still water depth, need to plot -phi + h/2
% where h is the uniform grid spacing
t = A(:,1);
eta = -A(:,2) + dx/2;
%shift = -dx/2.0; %max(B(:,2)) - max(y)
plot(TIME,y,'k-', t, eta,'r-');

% Legend
alpha = 2.0;
xlabel('Time');
ylabel('Elevation');
legend('Analytical Stokes 2nd', 'Numerical');
title(['x = ' num2str(x) ', alpha = ' num2str(alpha)]);