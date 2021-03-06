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

N = 50;

x = 3.8;

W = importdata('irregular_wave.txt');
Amp = W(:,1);
Ohm = W(:,2);
K  = W(:,3);
Ph = W(:,4);

D = importdata('probe_0');
A = D.data;
TIME = A(:,1);


y = zeros(1, length(TIME));
for j = 1: length(TIME)
    t = TIME(j);
    for i = 1:N
        y(j) = y(j) + Amp(i)*cos(K(i)*x - Ohm(i)*t + Ph(i));
    end
end


% If comparing from a cell above the still water depth, need to plot -phi + h/2
% where h is the uniform grid spacing
t = A(:,1);
eta = -A(:,2) + A(1,2);
plot(TIME,y,'k-', t, eta,'r-','linewidth',3);

% Legend
xlabel('Time');
ylabel('Elevation');
set(gca,'FontName','arial','fontSize',25,'fontWeight','normal','position',[0.165 0.19 0.8 0.75])
hleg = legend('Analytical irregular wave', 'IBAMR');
legend boxoff
set(hleg,'Position',[0.45, 0.2, 0.5, 0.1],'FontName','arial','fontSize',20);
%title(['x = ' num2str(x)  ' for Probe = ' num2str(0)]);

