clear all; close all; clc;

horizon = [3.015 2.015 1.015];


a = [17 33 65 97];
b = [9 17 33 49];

nstep = [200000 400000 136000 142000; ...
        200000 400000 387000 132000; ...
        200000 400000 270000 1200000];
ntime = [500 500 85 59.1667; ...
        500 500 241.875 55; ...
        500 500 168.75 500];

for l = 1:3
for k = 1:4

%     nstep(k) = 60000*(b(k)-1)/4;
%     ntime(k) = 300;
    clear A;
    A = load(['P=0.49995/d=' num2str(horizon(l)) '/' num2str((a(k)-1)/10) '/data/D_' num2str(nstep(l,k)) '_' num2str(ntime(l,k))]);

    length = a(k);
    num_rows = b(k);
    num_depth = 1;

    X = zeros(length, num_rows);
    Y = zeros(length, num_rows);
    U = zeros(length, num_rows);
    V = zeros(length, num_rows);
    D = zeros(length, num_rows);

    for j = 1:num_rows
        for i = 1:length
            idx = (j-1)*length*num_depth + i;
            X(i ,j) = A(idx, 1);
            Y(i ,j) = A(idx, 2);
            U(i ,j) = A(idx, 3);
            V(i ,j) = A(idx, 4);
            D(i ,j) = A(idx, 5);
            Xt(i,j) = A(idx, 1) + A(idx, 3);
            Yt(i,j) = A(idx, 2) + A(idx, 4);
        end
    end

    xpi = [X(1,:) X(2:length-1,num_rows)' X(length,num_rows:-1:1) X(length-1:-1:2,1)'];
    ypi = [Y(1,:) Y(2:length-1,num_rows)' Y(length,num_rows:-1:1) Y(length-1:-1:2,1)'];
    ini_area = polyarea(xpi,ypi);

    xp = [Xt(1,:) (Xt(2:length-1,num_rows))' Xt(length,num_rows:-1:1) (Xt(length-1:-1:2,1))'];
    yp = [Yt(1,:) (Yt(2:length-1,num_rows))' Yt(length,num_rows:-1:1) (Yt(length-1:-1:2,1))'];
    area(l,k) = (polyarea(xp,yp)-ini_area)/ini_area;

    total_nodes(l,k) = length * num_rows;
    y_disp(l,k) = A(length*num_rows - (length-1)/2,4);
end


blue    = [0.25,0.3,0.65];          % blue color : [0.25,0.3,0.65]
figure(2*l-1)
ppp(1) = plot(total_nodes(l,:),y_disp(l,:),'-o','color',blue                         ...
            , 'linewidth'       , 2.5                           ...
            , 'MarkerSize'      , 7                             ...
            , 'MarkerFaceColor' , blue                         ...
            , 'MarkerEdgeColor' , blue                         ...
            , 'MarkerIndices'   , 1:1:4 );
ylim([-4.5 0])
xlim([0 4900])
xlabel('DoF' , 'Interpreter', 'latex'   ...
      , 'FontName'       , 'arial'                  ...
      , 'fontSize'       , 25                       ...
      , 'fontWeight'     , 'bold')
ylabel('y-displacement (cm)' , 'Interpreter', 'latex'   ...
      , 'FontName'       , 'arial'                  ...
      , 'fontSize'       , 25                       ...
      , 'fontWeight'     , 'bold')
set(gca, 'FontName'      , 'arial'                  ...
       , 'fontSize'      , 18                       ...
       , 'fontWeight'    , 'normal'                 ...
       , 'position'      , [0.165 0.19 0.8 0.75]);
hleg = legend(ppp(1)  ,{'$\delta = 3.015 \Delta x$'}                ...
       , 'Interpreter', 'latex'                   ...
       , 'FontName'   , 'arial'                   ...
       , 'fontSize'   , 15                        ...
       , 'Location'   , 'northeast');
legend boxoff
title('\nu = 0.49995')

h = figure(2*l-1);
fig_props.figW = 16;
fig_props.figH = 12;
set(h, 'units'        , 'centimeters'     ...
     , 'color'        , [1 1 1]           ...
     , 'position'     , [1 1 fig_props.figW fig_props.figH] ...
     , 'PaperUnits'   , 'centimeters'     ...
     , 'PaperSize'    , [fig_props.figW fig_props.figH]     ...
     , 'PaperPosition', [0 0 fig_props.figW fig_props.figH])
print(h,'-dpdf','plot_name');

blue    = [0.25,0.3,0.65];          % blue color : [0.25,0.3,0.65]
figure(2*l)
ppp(1) = plot(total_nodes(l,:),area(l,:),'-o','color',blue                         ...
            , 'linewidth'       , 2.5                           ...
            , 'MarkerSize'      , 7                             ...
            , 'MarkerFaceColor' , blue                         ...
            , 'MarkerEdgeColor' , blue                         ...
            , 'MarkerIndices'   , 1:1:4 );
ylim([-0.1 0.1])
xlim([0 4900])
xlabel('DoF' , 'Interpreter', 'latex'   ...
      , 'FontName'       , 'arial'                  ...
      , 'fontSize'       , 25                       ...
      , 'fontWeight'     , 'bold')
ylabel('Volume Change (%)' , 'Interpreter', 'latex'   ...
      , 'FontName'       , 'arial'                  ...
      , 'fontSize'       , 25                       ...
      , 'fontWeight'     , 'bold')
set(gca, 'FontName'      , 'arial'                  ...
       , 'fontSize'      , 18                       ...
       , 'fontWeight'    , 'normal'                 ...
       , 'position'      , [0.165 0.19 0.8 0.75]);
hleg = legend(ppp(1)  ,{'$\delta = 3.015 \Delta x$'}                ...
       , 'Interpreter', 'latex'                   ...
       , 'FontName'   , 'arial'                   ...
       , 'fontSize'   , 15                        ...
       , 'Location'   , 'northeast');
legend boxoff
title('\nu = 0.49995')

h = figure(2*l);
fig_props.figW = 16;
fig_props.figH = 12;
set(h, 'units'        , 'centimeters'     ...
     , 'color'        , [1 1 1]           ...
     , 'position'     , [1 1 fig_props.figW fig_props.figH] ...
     , 'PaperUnits'   , 'centimeters'     ...
     , 'PaperSize'    , [fig_props.figW fig_props.figH]     ...
     , 'PaperPosition', [0 0 fig_props.figW fig_props.figH])
print(h,'-dpdf','plot_name');
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

% xpi = [X(1,:) X(2:length-1,num_rows)' X(length,num_rows:-1:1) X(length-1:-1:2,1)'];
% ypi = [Y(1,:) Y(2:length-1,num_rows)' Y(length,num_rows:-1:1) Y(length-1:-1:2,1)'];
% ini_area = polyarea(xpi,ypi);

% for i = 1:10
%     nstep(i) = i*1000;
%     ntime(i) = i*1.25;
%     B = load(['data/D_' num2str(nstep(i)) '_' num2str(ntime(i))]);
% 
%     for j = 1:num_rows
%         for k = 1:length
%             idx = (j-1)*length*num_depth + k;
%             Xt(k ,j) = B(idx, 1) + B(idx, 3);
%             Yt(k ,j) = B(idx, 2) + B(idx, 4);
%         end
%     end
%     xp = [Xt(1,:) (Xt(2:length-1,num_rows))' Xt(length,num_rows:-1:1) (Xt(length-1:-1:2,1))'];
%     yp = [Yt(1,:) (Yt(2:length-1,num_rows))' Yt(length,num_rows:-1:1) (Yt(length-1:-1:2,1))'];
%     area(i) = (polyarea(xp,yp)-ini_area)/ini_area;
% 
%     y_disp(i) = B(length*num_rows - (length-1)/2,4);
% end
% 
% figure(4)
% plot([0 ntime],[0 y_disp],'LineWidth',3.0)
% xlabel('time')
% ylabel('U_y')
% 
% figure(5)
% plot([0 ntime],[0 area],'LineWidth',3.0)
% xlabel('time')
% ylabel('A')

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