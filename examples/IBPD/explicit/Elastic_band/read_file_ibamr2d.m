clear all; close all; clc;

length = 5;
num_rows = 40;
num_depth = 1;

B = load('data/D_1000_0.025');

    for k = 1:length
        for j = 1:num_rows
            idx = (j-1)*length + k;
            X(j,k) = B(idx, 1);
            Y(j,k) = B(idx, 2);
        end
    end

yyy1 = 0.1 * ones(1,length);
yyy2 = 0.9 * ones(1,length);
    
xpi = [X(length,:) X(length:num_rows-length+1,length)' X(num_rows-length+1,length:-1:1) X(num_rows-length+1:-1:length,1)'];
ypi = [yyy1 Y(length:num_rows-length+1,length)' yyy2 Y(num_rows-length+1:-1:length,1)'];
ini_area = polyarea(xpi,ypi);

figure(3)
scatter(xpi,ypi)

for i = 290
    nstep(i) = i*1000;
    ntime(i) = i*0.025;
    B = load(['data/D_' num2str(nstep(i)) '_' num2str(ntime(i))]);

    for k = 1:length
        for j = 1:num_rows
            idx = (j-1)*length + k;
            Xt(j,k) = B(idx, 1) + B(idx,3);
            Yt(j,k) = B(idx, 2) + B(idx,4);
            U(j,k) = B(idx, 3);
            V(j,k) = B(idx, 4);
        end
    end
    
    xp = [ Xt(length,1) X(length,:) Xt(length,length) Xt(length:num_rows-length+1,length)' Xt(num_rows-length+1,length) X(num_rows-length+1,length:-1:1) Xt(num_rows-length+1,1) Xt(num_rows-length+1:-1:length,1)'];
    yp = [Yt(length,1) yyy1 Yt(length,length) Yt(length:num_rows-length+1,length)' Yt(num_rows-length+1,length) yyy2 Yt(num_rows-length+1,1) Yt(num_rows-length+1:-1:length,1)'];
    area(i) = (polyarea(xp,yp)-ini_area)/ini_area * 100;
    
    figure(2)
    scatter(xp,yp);
    
    x_disp(i) = (U(num_rows/2,length) + U(num_rows/2 + 1,length))/2;
end

figure(4)
plot([0 ntime],[0 area],'LineWidth',3.0)
xlabel('time')
ylabel('Volume change (%)')

figure(5)
plot([0 ntime],[0 x_disp],'LineWidth',3.0)
xlabel('time')
ylabel('U_x')