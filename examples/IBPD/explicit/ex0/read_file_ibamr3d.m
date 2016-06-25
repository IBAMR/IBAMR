% Filename : read_file_ibamr3d.m % Created on 24 Jun 2016 by Amneet Bhalla %
    All rights reserved.

    clear all;
clc;
cla;

A = load('data/D_7000_0.0007');

length = 151;
num_rows = 76;
num_depth = 3;
X = zeros(num_rows, length);
Y = zeros(num_rows, length);
U = zeros(num_rows, length);
V = zeros(num_rows, length);
D = zeros(num_rows, length);

for i = 1:length
    for j = 1:num_rows
        idx = (i-1)*num_rows*num_depth + num_rows + j;
X(j, i) = A(idx, 1);
Y(j, i) = A(idx, 2);
U(j, i) = A(idx, 3);
V(j, i) = A(idx, 4);
D(j, i) = A(idx, 5);
end end

    figure(1);
pcolor(X, Y, U);
title('u_x');
shading flat colorbar

    figure(2);
pcolor(X, Y, V);
title('u_y');
shading flat colorbar

    figure(5) plot(X(3,
                     :),
                   U(3,
                     :),
                   'k*-',
                   X(3,
                     :),
                   (200E6 / 200E9) * (X(3,
                                        :) -
                                      0.5),
                   '--');
xlabel('x') ylabel('u_x') legend('PD solution', 'Analytical');
hold on;

figure(3) plot(Y( :, 75), V( :, 75), 'k*-', Y( :, 75), -(1.0 / 3.0) * (200E6 / 200E9) * (Y( :, 75) - 0.25), '--');
xlabel('y') ylabel('u_y') legend('PD solution', 'Analytical');
hold on;
% %

    A = load('data/D_9500_0.00095');
X = zeros(num_rows, length);
Y = zeros(num_rows, length);
U = zeros(num_rows, length);
V = zeros(num_rows, length);
D = zeros(num_rows, length);

for i = 1:length
    for j = 1:num_rows
        idx = (i-1)*num_rows*num_depth + num_rows + j;
X(j, i) = A(idx, 1);
Y(j, i) = A(idx, 2);
U(j, i) = A(idx, 3);
V(j, i) = A(idx, 4);
D(j, i) = A(idx, 5);
end end

    figure(3);
plot(Y( :, 75), V( :, 75), 'g*-', Y( :, 75), -(1.0 / 3.0) * (200E6 / 200E9) * (Y( :, 75) - 0.25), '--');
hold off;

figure(4);
pcolor(X, Y, U);
shading flat colorbar

    figure(5) plot(X(3,
                     :),
                   U(3,
                     :),
                   'c*-',
                   X(3,
                     :),
                   (200E6 / 200E9) * (X(3,
                                        :) -
                                      0.5),
                   '--');
xlabel('x') ylabel('u_x') legend('PD solution', 'Analytical');
hold off
