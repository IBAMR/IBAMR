clear all; close all; clc;

A = load('data/D_50000_5');

length = 13;
num_rows = 3;
num_depth = 3;

total_nodes(1) = length * num_rows * num_depth;

Top_center_pt_idx = (length-1) * num_rows * num_depth + (num_rows * num_depth + 1)/2;
Max_disp(1) = A(Top_center_pt_idx,4);

B = load('data/D_100000_5');

length = 25;
num_rows = 5;
num_depth = 5;

total_nodes(2) = length * num_rows * num_depth;

Top_center_pt_idx = (length-1) * num_rows * num_depth + (num_rows * num_depth + 1)/2;
Max_disp(2) = B(Top_center_pt_idx,4);

C = load('data/D_200000_5');

length = 49;
num_rows = 9;
num_depth = 9;

total_nodes(3) = length * num_rows * num_depth;

Top_center_pt_idx = (length-1) * num_rows * num_depth + (num_rows * num_depth + 1)/2;
Max_disp(3) = C(Top_center_pt_idx,4);


blue    = [0.25,0.3,0.65];          % blue color : [0.25,0.3,0.65]
figure(2)
ppp(1) = plot(total_nodes,Max_disp,'-o','color',blue                         ...
            , 'linewidth'       , 1.5                           ...
            , 'MarkerSize'      , 5                             ...
            , 'MarkerFaceColor' , blue                         ...
            , 'MarkerEdgeColor' , blue                         ...
            , 'MarkerIndices'   , 1:1:3 );
ylim([-1.5 0.9])
xlim([0 8000])
xlabel('DoF' , 'Interpreter', 'latex'   ...
      , 'FontName'       , 'arial'                  ...
      , 'fontSize'       , 25                       ...
      , 'fontWeight'     , 'bold')
ylabel('Disp.(cm)' , 'Interpreter', 'latex'   ...
      , 'FontName'       , 'arial'                  ...
      , 'fontSize'       , 25                       ...
      , 'fontWeight'     , 'bold')
set(gca, 'FontName'      , 'arial'                  ...
       , 'fontSize'      , 18                       ...
       , 'fontWeight'    , 'normal'                 ...
       , 'position'      , [0.165 0.19 0.8 0.75]);
hleg = legend(ppp(1)  ,{'IPD'}                ...
       , 'Interpreter', 'latex'                   ...
       , 'FontName'   , 'arial'                   ...
       , 'fontSize'   , 15                        ...
       , 'Location'   , 'northeast');
legend boxoff

h = figure(2);
fig_props.figW = 20;
fig_props.figH = 16;
set(h, 'units'        , 'centimeters'     ...
     , 'color'        , [1 1 1]           ...
     , 'position'     , [1 1 fig_props.figW fig_props.figH] ...
     , 'PaperUnits'   , 'centimeters'     ...
     , 'PaperSize'    , [fig_props.figW fig_props.figH]     ...
     , 'PaperPosition', [0 0 fig_props.figW fig_props.figH])
print(h,'-dpdf','plot_name');

% X = zeros(num_rows, num_depth, length);
% Y = zeros(num_rows, num_depth, length);
% Z = zeros(num_rows, num_depth, length);
% U = zeros(num_rows, num_depth, length);
% V = zeros(num_rows, num_depth, length);
% W = zeros(num_rows, num_depth, length);
% D = zeros(num_rows, num_depth, length);
% 
% for k = 1:length
%     for i = 1:num_depth
%         for j = 1:num_rows
%             idx = (k-1)*num_rows*num_depth + (i-1)*num_rows + j;
%             X(i,j,k) = A(idx, 1);
%             Y(i,j,k) = A(idx, 2);
%             Z(i,j,k) = A(idx, 3);
%             U(i,j,k) = A(idx, 4);
%         end
%     end
% end

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
