clear all; close all; clc;

A = load('data/D_5000_0.005');

length = 101;
num_rows = 51;
num_depth = 1;

X = zeros(num_rows, length);
Y = zeros(num_rows, length);
U = zeros(num_rows, length);
V = zeros(num_rows, length);
D = zeros(num_rows, length);

for i = 1:length
    for j = 1:num_rows
        idx = (i-1)*num_rows*num_depth + j;
        X(j, i) = A(idx, 1);
        Y(j, i) = A(idx, 2);
        U(j, i) = A(idx, 3);
        V(j, i) = A(idx, 4);
        D(j, i) = A(idx, 5);
        p(j+num_rows*(i-1),1) = A(idx, 1) + A(idx, 3);
        p(j+num_rows*(i-1),2) = A(idx, 2) + A(idx, 4);
    end
end

figure(1);
pcolor(X, Y, U);
title('u_x');
shading flat
colorbar

figure(2);
pcolor(X, Y, V);
title('u_y');
shading flat
colorbar

figure(3)
scatter(p(:,1),p(:,2),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on

% n1 = (num_rows+1)/2;
% n2 = (length+1)/2;
% 
% figure(4)
% plot(X(n1,:),U(n1,:),'k*-', X(n1,:),(1E3 / 1E6) * (X(n1,:) -0.5),'--');
% xlabel('x')
% ylabel('u_x')
% legend('PD solution', 'Analytical');
% hold on;
% 
% figure(5)
% plot(Y( :, n2), V( :, n2), 'k*-', Y( :, n2), -0.4 * (1E3 / 1E6) * (Y(:,n2) - 0.25), '--');
% xlabel('y')
% ylabel('u_y')
% legend('PD solution', 'Analytical');
% hold on
