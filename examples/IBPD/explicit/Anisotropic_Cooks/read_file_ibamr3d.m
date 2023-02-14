% Filename : read_file_ibamr3d.m % Created on 24 Jun 2016 by Amneet Bhalla %
% All rights reserved.

clear all; close all; clc;

A = load('data/D_6000_0.006');

length = 101;
num_rows = 51;
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

n1 = (num_rows+1)/2;
n2 = (length+1)/2;

figure(4)
plot(X(n1,:),U(n1,:),'k*-', X(n1,:),(1E3 / 1E6) * (X(n1,:) -0.5),'--');
xlabel('x')
ylabel('u_x')
legend('PD solution', 'Analytical');
hold on;

figure(5)
plot(Y( :, n2), V( :, n2), 'k*-', Y( :, n2), -0.4 * (1E3 / 1E6) * (Y(:,n2) - 0.25), '--');
xlabel('y')
ylabel('u_y')
legend('PD solution', 'Analytical');
hold on


% % %Basic plot colors.
% % black   = [0 0 0];                  % black color: [0 0 0]
% % red     = [0.95,0.1,0.15];          % red  color : [0.95,0.1,0.15]
% % green   = [0.5,0.8,0.35];           % green color: [0.5,0.8,0.35]
% % blue    = [0.25,0.3,0.65];          % blue color : [0.25,0.3,0.65]
% % mustard = [1.0, 0.8, 0.4];          % mustard color: [1.0, 0.8, 0.4]
% % magenta = [0.8500 0.3250 0.0980];   % magenta color: [0.8500 0.3250 0.0980]
% %
% % figure(5)
% % ppp(1) = plot(X(n1,:),U(n1,:),'-', 'color', red                         ...
% %             , 'linewidth'       , 1.5                           ...
% %             , 'MarkerSize'      , 5                             ...
% %             , 'MarkerFaceColor' , red                         ...
% %             , 'MarkerEdgeColor' , red                         ...
% %             , 'MarkerIndices'   , 1:1:length );
% % hold on;
% % ppp(2) = plot(X(n1,:),(1E4 / 1E6) * (X(n1,:) -0.5),'--','color',blue                         ...
% %             , 'linewidth'       , 1.5                           ...
% %             , 'MarkerSize'      , 5                             ...
% %             , 'MarkerFaceColor' , blue                         ...
% %             , 'MarkerEdgeColor' , blue                         ...
% %             , 'MarkerIndices'   , 1:1:length );
% % hold off;
% % xlabel('$\mathbf{x}$' , 'Interpreter', 'latex'   ...
% %       , 'FontName'       , 'arial'                  ...
% %       , 'fontSize'       , 25                       ...
% %       , 'fontWeight'     , 'bold')
% % ylabel('$\mathbf{U_x}$' , 'Interpreter', 'latex'   ...
% %       , 'FontName'       , 'arial'                  ...
% %       , 'fontSize'       , 25                       ...
% %       , 'fontWeight'     , 'bold')
% % set(gca, 'FontName'      , 'arial'                  ...
% %        , 'fontSize'      , 18                       ...
% %        , 'fontWeight'    , 'normal'                 ...
% %        , 'position'      , [0.165 0.19 0.8 0.75]);
% % hleg = legend(ppp(1:2)  ,{'Peridynamic solution', 'Analytic solution'}                ...
% %        , 'Interpreter', 'latex'                   ...
% %        , 'FontName'   , 'arial'                   ...
% %        , 'fontSize'   , 15                        ...
% %        , 'Location'   , 'northwest');
% % legend boxoff
% %
% % h = figure(5);
% % fig_props.figW = 24;
% % fig_props.figH = 16;
% % set(h, 'units'        , 'centimeters'     ...
% %      , 'color'        , [1 1 1]           ...
% %      , 'position'     , [1 1 fig_props.figW fig_props.figH] ...
% %      , 'PaperUnits'   , 'centimeters'     ...
% %      , 'PaperSize'    , [fig_props.figW fig_props.figH]     ...
% %      , 'PaperPosition', [0 0 fig_props.figW fig_props.figH])
% % print(h,'-dpdf','plot_name');
