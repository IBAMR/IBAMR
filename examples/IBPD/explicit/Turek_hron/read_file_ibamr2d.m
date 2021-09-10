clear all; close all; clc;

length = 125;
num_rows = 7;
num_depth = 1;


for i = 1:30
    nstep(i) = i*1000;
    ntime(i) = i*0.003;
    B = load(['data/D_' num2str(nstep(i)) '_' num2str(ntime(i))]);

    for k = 1:length
        for j = 1:num_rows
            idx = (j-1)*length + k;
            X(j,k) = B(idx, 1);
            Y(j,k) = B(idx, 2);
            U(j,k) = B(idx, 3);
            V(j,k) = B(idx, 4);
        end
    end
    
    y_disp(i) = V((num_rows+1)/2,length);
    x_disp(i) = U((num_rows+1)/2,length);
end

figure(4)
plot([0 ntime],[0 y_disp],'LineWidth',3.0)
xlabel('time')
ylabel('U_y')

figure(5)
plot([0 ntime],[0 x_disp],'LineWidth',3.0)
xlabel('time')
ylabel('U_x')
