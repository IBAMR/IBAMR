clear all; close all; clc;

A = load('data/D_1000_0.2');

length = 49;
num_rows = 45;
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

% figure(1);
% pcolor(X, Y, U);
% title('u_x');
% shading flat
% colorbar
% 
% figure(2);
% pcolor(X, Y, V);
% title('u_y');
% shading flat
% colorbar
% 
% figure(3)
% scatter(p(:,1),p(:,2),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% grid on
% 
% n1 = (num_rows)/2;
% n2 = (length+1)/2;
% 
% figure(4)
% plot(X(n1,:),U(n1,:),'k*-');
% xlabel('x')
% ylabel('u_x')
% legend('PD solution');
% hold on;
% 
% figure(5)
% plot(Y( :, n2), V( :, n2), 'k*-');
% xlabel('y')
% ylabel('u_y')
% legend('PD solution');
% hold on

xpi = [X(1,:) X(2:num_rows-1,length)' X(num_rows,length:-1:1) X(num_rows-1:-1:2,1)'];
ypi = [Y(1,:) Y(2:num_rows-1,length)' Y(num_rows,length:-1:1) Y(num_rows-1:-1:2,1)'];
ini_area = polyarea(xpi,ypi);

for i = 1:100
    nstep(i) = i*1000;
    ntime(i) = i*0.2;
    B = load(['data/D_' num2str(nstep(i)) '_' num2str(ntime(i))]);

    for j = 1:length
        for k = 1:num_rows
            idx = (j-1)*num_rows*num_depth + k;
            Xt(k ,j) = B(idx, 1) + B(idx, 3);
            Yt(k ,j) = B(idx, 2) + B(idx, 4);
        end
    end
    
    xp = [Xt(1,:) Xt(2:num_rows-1,length)' Xt(num_rows,length:-1:1) Xt(num_rows-1:-1:2,1)'];
    yp = [Yt(1,:) Yt(2:num_rows-1,length)' Yt(num_rows,length:-1:1) Yt(num_rows-1:-1:2,1)'];
    area(i) = (polyarea(xp,yp)-ini_area)/ini_area;
    
    y_disp(i) = B(length*num_rows,4);
end

figure(4)
plot([0 ntime],[0 y_disp],'LineWidth',3.0)
xlabel('time')
ylabel('U_y')

figure(5)
plot([0 ntime],[0 area],'LineWidth',3.0)
xlabel('time')
ylabel('A')

% % n1 = (num_rows+1)/2;
% % n2 = (length+1)/2;
% % 
% % figure(4)
% % plot(X(n1,:),U(n1,:),'k*-', X(n1,:),(1E3 / 1E6) * (X(n1,:) -0.5),'--');
% % xlabel('x')
% % ylabel('u_x')
% % legend('PD solution', 'Analytical');
% % hold on;
% % 
% % figure(5)
% % plot(Y( :, n2), V( :, n2), 'k*-', Y( :, n2), -0.4 * (1E3 / 1E6) * (Y(:,n2) - 0.25), '--');
% % xlabel('y')
% % ylabel('u_y')
% % legend('PD solution', 'Analytical');
% % hold on
