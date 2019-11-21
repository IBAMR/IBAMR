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


DEPTH       =  2.35;

AMPLITUDE   =  0.0791044;

WAVELENGTH  =  2.93483;               

WAVENUMBER  =  2*pi/WAVELENGTH;

TIME_PERIOD_STOKES = 1.3711;

OMEGA_STOKES = 2*pi/TIME_PERIOD_STOKES;

GRAVITY     = 9.81;

TIME = linspace(0,40,1000);

WAVE_SPACE = [DEPTH, 2*AMPLITUDE]/(GRAVITY * TIME_PERIOD_STOKES^2) 

x     = 7.34375;

D = importdata('probe_0');
A = D.data;


% Calculate elevation by 5th order Stokes
[C,eta] = stokes_coefs(WAVENUMBER, DEPTH, AMPLITUDE, GRAVITY);

k = WAVENUMBER;
ka = k * AMPLITUDE;
g = GRAVITY;
OMEGA_5th = sqrt(g * k) * (C(1) + ka^2 * C(3) + ka^4 * C(5));
phase = k * x - OMEGA_5th * TIME;

y_5th = (eta(1) * cos(phase) + eta(2) * cos(2 * phase) + eta(3) * cos(3 * phase) + eta(4) * cos(4 * phase) + ...
            eta(5) * cos(5 * phase))/ k;

% Calculate elevation by 2nd order Stokes.
H = 2*AMPLITUDE;
y_2nd = H/2*cos(WAVENUMBER*x - OMEGA_STOKES*TIME) + ... 
    pi*H^2/(8*WAVELENGTH) * cosh(WAVENUMBER*DEPTH) * ...
    (2+cosh(2*WAVENUMBER*DEPTH))/(sinh(WAVENUMBER*DEPTH)^3) * ... 
    cos(2*(WAVENUMBER*x - OMEGA_STOKES*TIME));



% If comparing from a cell above the still water depth, need to plot -phi + h/2
% where h is the uniform grid spacing
t = A(:,1);
eta = -A(:,2) + A(1,2);
figure;
plot(TIME,y_5th,'k-', t, eta,'r-');


% Legend
alpha = 2.0;
xlabel('Time');
ylabel('Elevation');
legend('Analytical Stokes 5th', 'Analytical Stokes 2nd', 'Numerical');
title(['x = ' num2str(x) ', alpha = ' num2str(alpha)]);


function [C, eta] = stokes_coefs(k,d,a,g)

ka = k*a;
kd = k*d;
ka2 = ka^2;
ka3 = ka^3;
ka4 = ka^4;
ka5 = ka^5;

A = zeros(5,5);
B = zeros(5,5);
p = zeros(5,1);
C = zeros(5,1);
eta = zeros(5,1);

S = 1 / cosh(2 * kd);

% Terms to calculate velocity.
A(1,1) = 1 / sinh(kd);
A(2,2) = 3 * S * S / (2 * (1 - S) * (1 - S));
A(3,1) = (-4 - 20 * S + 10 * S * S - 13 * S * S * S) / (8 * sinh(kd) * pow((1 - S), 3));
A(3,3) = (-2 * S * S + 11 * S * S * S) / (8 * sinh(kd) * pow((1 - S), 3));
A(4,2) = (12 * S - 14 * S * S - 264 * S * S * S - 45 * S * S * S * S - 13 * S * S * S * S * S) / ...
            (24 * pow((1 - S), 5));
A(4,4) = (10 * S * S * S - 174 * S * S * S * S + 291 * pow(S, 5) + 278 * pow(S, 6)) / ...
            (48 * (3 + 2 * S) * pow((1 - S), 5));
A(5,1) = (-1184 + 32 * S + 13232 * S * S + 21712 * S * S * S + 20940 * S * S * S * S + 12554 * pow(S, 5) - ...
             500 * pow(S, 6) - 3341 * pow(S, 7) - 670 * pow(S, 8)) / ...
            (64 * sinh(kd) * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));
A(5,3) = (4 * S + 105 * S * S + 198 * S * S * S - 1376 * S * S * S * S - 1302 * pow(S, 5) - 117 * pow(S, 6) + ...
             58 * pow(S, 7)) / ...
            (32 * sinh(kd) * (3 + 2 * S) * pow((1 - S), 6));
A(5,5) = (-6 * S * S * S + 272 * S * S * S * S - 1552 * pow(S, 5) + 852 * pow(S, 6) + 2029 * pow(S, 7) + ...
             430 * pow(S, 8)) / ...
            (64 * sinh(kd) * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));

% Terms to calculate eta.
B(2,2) = 1 / tanh(kd) * (1 + 2 * S) / (2 * (1 - S));
B(3,1) = -3 * (1 + 3 * S + 3 * S * S + 2 * S * S * S) / (8 * pow((1 - S), 3));
B(4,2) = 1 / tanh(kd) * (6 - 26 * S - 182 * S * S - 204 * S * S * S - 25 * pow(S, 4) + 26 * pow(S, 5)) / ...
            (6 * (3 + 2 * S) * pow((1 - S), 4));
B(4,4) = 1 / tanh(kd) * (24 + 92 * S + 122 * S * S + 66 * S * S * S + 67 * pow(S, 4) + 34 * pow(S, 5)) / ...
            (24 * (3 + 2 * S) * pow((1 - S), 4));
B(5,3) = 9 *(132 + 17 * S - 2216 * S * S - 5897 * S * S * S - 6292 * pow(S, 4) - 2687 * pow(S, 5) + ...
             194 * pow(S, 6) + 467 * pow(S, 7) + 82 * pow(S, 8)) / ...
            (128 * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));
B(5,5) = 5 *(300 + 1579 * S + 3176 * S * S + 2949 * S * S * S + 1188 * pow(S, 4) + 675 * pow(S, 5) + ...
             1326 * pow(S, 6) + 827 * pow(S, 7) + 130 * pow(S, 8)) / ...
            (384 * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));

% Wave dispersion
C(1) = sqrt(tanh(kd));
C(3) = C(1) * (2 + 7 * S * S) / (4 * pow((1 - S), 2));
C(5) = C(1) * (4 + 32 * S - 116 * S * S - 400 * S * S * S - 71 * pow(S, 4) + 146 * pow(S, 5)) / ...
         (32 * pow((1 - S), 5));

% Velocity coefficients.
alpha = C(1) * sqrt(g / pow(k, 3));
p(1) = alpha * (ka * A(1,1) + ka3 *  A(3,1) + ka5 * A(5,1));
p(2) = alpha * (ka2 * A(2,2) + ka4 * A(4,2));
p(3) = alpha * (ka3 * A(3,3) + ka5 * A(5,3));
p(4) = alpha * (ka4 * A(4,4));
p(5) = alpha * (ka5 * A(5,5));

% Component wave eta.
eta(1) = ka + ka3 * B(3,1) - ka5 * (B(5,3) + B(5,5));
eta(2) = ka2 * B(2,2) + ka4 * B(4,2);
eta(3) = -ka3 * B(3,1) + ka5 * B(5,3);
eta(4) = ka4 * B(4,4);
eta(5) = ka5 * B(5,5);

end

function c = pow(a,b)

c = a^b; 

end

