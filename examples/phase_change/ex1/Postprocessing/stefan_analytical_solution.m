% Solving transcendental equation and plotting the analytical interface
% position as a function of time for one-dimensional Stefan problem with
% expansion (\rho_s/\rho_l < 1).

clc;
clear all;
close all;


%Basic plot colors.
black   = [0 0 0];                  % black color: [0 0 0]
red     = [0.95,0.1,0.15];          % red  color : [0.95,0.1,0.15]
green   = [0,0.6,0];                % green color: [0.0,0.6,0.0]
blue    = [0.25,0.3,0.65];          % blue color : [0.25,0.3,0.65]
mustard = [1.0, 0.8, 0.4];          % mustard color: [1.0, 0.8, 0.4]
magenta = [0.8500 0.3250 0.0980];   % magenta color: [0.8500 0.3250 0.0980]
orange  = [1, 0.5, 0];
yellow = [1, 1, 0];

% Properties of liquid and solid phases
cp_s = 910;
k_s  = 211.0;
k_l  = 91.0;
cp_l = 1042.4;
L    = 3.8384e5;
rho_l  = 2700;
rho_s  = 500;
T_melt = 933.6;

% Left and right bcs
T_liquid = 973.6;
T_0      = 298.6;

% simulation time
dt = 1e-3;
end_time = 10.0;
counter = 1.0;
time = dt:dt:end_time;

alpha_s = k_s/(rho_s*cp_s);
alpha_l = k_l/(rho_l*cp_l);

dT_s = T_melt - T_0;
dT_l = T_melt - T_liquid;
R_rho = rho_s/rho_l;
diffusivity_ratio = alpha_l/alpha_s;
alpha_st = sqrt(alpha_s);
alpha_lt = sqrt(alpha_l);
for i=1:length(time)
func = @(lambda) k_s*(dT_s*exp(-lambda^2*diffusivity_ratio)/ (sqrt(pi).*alpha_st*erf(lambda*sqrt(diffusivity_ratio)))) ...
    + k_l*(dT_l*exp(-(lambda*R_rho)^2)/ (sqrt(pi).*alpha_lt*erfc(lambda*R_rho))) ...
    - rho_s*lambda*sqrt(alpha_l)*(L - 0.5*(1-R_rho^2)*lambda^2*alpha_l/time(i));

%  Make sure the solution is physical. If not try some other initial
%  conditions.
lambda_ke(i) = fzero(func,2.0);
x_0 = 0.0;
x_analytical_ke(i) = x_0 + 2*lambda_ke(i)*sqrt(time(i)*alpha_l);
end

plot(time,x_analytical_ke,'linewidth',3)