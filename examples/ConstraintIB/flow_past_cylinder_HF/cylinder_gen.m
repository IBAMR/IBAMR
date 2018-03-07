%% Filename : Cylinder2d.m

clf
clear all;
clc;
%% Mesh parameters.
Lx = 18.0; Ly = 12.0;
Nx = 900*8; Ny = 600*8;
dx = Lx/Nx; dy = Ly/Ny;

R = 0.5;

%% disk parameters.
X_com = 0.0; Y_com = 0.0;
h = dx/sqrt(2);
% Radius = 0.5-h;
% num_pts_x = ceil(2*Radius/dx)+1;
% num_pts_y = ceil(2*Radius/dy)+1;

%% write out the points.
% idx = 0;
% for i = 1:num_pts_x
%     x = X_com + ( (i-1)*dx - Radius );
%     
%     for j = 1:num_pts_y
%         y = Y_com + ( (j-1)*dy - Radius );
%         
%         if( ((x- X_com)^2 + (y-Y_com)^2) <= Radius^2 )
%             idx = idx+1;
%             X_array(idx) = x; Y_array(idx) = y;
%         end
%     end
% end

idx = 0;
circ_idx_begin = idx; 
radius = R; 
while radius > 0
circumference = 2 * pi * radius; 
num_pts_surf = round(circumference/dx);
h = circumference / num_pts_surf;
for i = 1:(num_pts_surf+1*0)
    s = (i-1) * h; 
    idx = idx+1;
    X_array(idx) = radius * cos(s / radius); 
    Y_array(idx) = radius * sin(s / radius); 
end
plot(X_array, Y_array, '.');
axis('equal')
hold on 

if radius == R
    circ_idx_end = idx-1;
end

radius = radius - h; 
end
% if (Nx ~= 64*4)
%     idx = idx+1;
%     X_array(idx) = 0.0; 
%     Y_array(idx) = 0.0; 
% end

XCOM = sum(X_array)/length(X_array);
YCOM = sum(Y_array)/length(Y_array);

X_array = X_array - XCOM;
Y_array = Y_array - YCOM;

slave_idx = [(circ_idx_begin+1):circ_idx_end, circ_idx_begin];
master_idx = circ_idx_begin:circ_idx_end; 
elem = [slave_idx; master_idx];

%% plot the disk
clf
plot(X_array, Y_array, '.');
axis('equal')
hold on 
plot(X_array(circ_idx_begin+1), Y_array(circ_idx_begin+1), 'ro')
plot(X_array(circ_idx_end+1), Y_array(circ_idx_end+1), 'ro')
grid on 
set(gca,'xtick',-1:dx:1)
set(gca,'ytick',-1:dy:1)
axis('equal')

for i = 1:length(elem)
plot([X_array(elem(1,i)+1), X_array(elem(2,i)+1)], ...
     [Y_array(elem(1,i)+1), Y_array(elem(2,i)+1)], '-')
end

% lifted surface
X_mid = 0.5 * (X_array(elem(1,:)+1) + X_array(elem(2,:)+1));
Y_mid = 0.5 * (Y_array(elem(1,:)+1) + Y_array(elem(2,:)+1));
plot(X_mid, Y_mid, 'o', 'MarkerSize', 4);
t_x = X_array(elem(2,:)+1) - X_array(elem(1,:)+1);
t_y = Y_array(elem(2,:)+1) - Y_array(elem(1,:)+1);
ds = sqrt(t_x.^2 + t_y.^2);
% quiver(X_mid, Y_mid, t_x./ds, t_y./ds)
n_x = -t_y./ds;
n_y = t_x./ds;
quiver(X_mid, Y_mid, n_x, n_y)
X_lifted = X_mid + 2*dx*n_x;
Y_lifted = Y_mid + 2*dy*n_y;
plot(X_lifted, Y_lifted, 's')
xlabel('x')
ylabel('y')

%% write the coordinates in the file
fid = fopen('cylinder2d.vertex','wt');
fprintf(fid,'%d\n', length(X_array));
for i = 1:length(X_array)
    fprintf(fid,'%8.16f\t%8.16f\n',X_array(i),Y_array(i));
end
fclose(fid); 

%% write the elements for ghost surface in the file
fid = fopen('cylinder2d.elem','wt');
fprintf(fid,'%d\n', length(elem));

for i = 1:length(elem)
    fprintf(fid,'%d\t%d\n',elem(2,i),elem(1,i));
end

%% rotations

omega = -1; 
theta = 0.0; 
Xo = X_array;
Yo = Y_array; 
N = length(Xo); 
Xcm = sum(Xo) / N; 
Ycm = sum(Yo) / N; 
U = zeros(1,N);
V = zeros(1,N); 

dt = 1; 
clf
for n = 0:floor(pi/dt)
    t = n * dt; 
    theta = theta + dt * omega; 
    
    C = cos(theta); 
    S = sin(theta); 
    
    X = (Xo - Xcm) * C - (Yo - Ycm) * S;
    Y = (Xo - Xcm) * S + (Yo - Ycm) * C;
    
    U = -omega * Y;
    V = omega * X;
    
    plot(X, Y, '.')
    hold on 
    quiver(X, Y, U, V)
    hold off
    title(strcat(['t = ', num2str(t)]))
    grid minor
    axis equal
    pause(0.01)
end

%% ghost surfaces
for j = [1:3, 8:8:24]
radius = R + j*dx; 
j*dx
circumference = 2 * pi * radius; 
num_pts_surf = round(circumference/dx);
h = circumference / num_pts_surf;
X_ghost = zeros(num_pts_surf, 1); 
Y_ghost = zeros(num_pts_surf, 1); 
for i = 1:(num_pts_surf+1*0)
    s = (i-1) * h; 
    X_ghost(i) = radius * cos(s / radius); 
    Y_ghost(i) = radius * sin(s / radius); 
end
plot(X_ghost, Y_ghost, '.');
axis('equal')
hold on 

%% write the coordinates in the file
fid = fopen(strcat(['cylinder2d_', num2str(j), 'h.ghost']),'wt');
fprintf(fid,'%d\n', num_pts_surf);
for i = 1:num_pts_surf
    fprintf(fid,'%8.16f\t%8.16f\n',X_ghost(i),Y_ghost(i));
end
fclose(fid); 

%% write the spring elements in the file

slave_idx = [1:(num_pts_surf-1), 0];
master_idx = 0:(num_pts_surf-1); 
elem = [slave_idx; master_idx];

fid = fopen(strcat(['cylinder2d_', num2str(j), 'h.elem']),'wt');
fprintf(fid,'%d\n', length(master_idx));
for i = 1:length(master_idx)
    fprintf(fid,'%d\t%d\n', master_idx(i), slave_idx(i));
end
fclose(fid); 

% lifted surface
X_mid = 0.5 * (X_ghost(elem(1,:)+1) + X_ghost(elem(2,:)+1));
Y_mid = 0.5 * (Y_ghost(elem(1,:)+1) + Y_ghost(elem(2,:)+1));
plot(X_mid, Y_mid, 'o', 'MarkerSize', 4);
t_x = X_ghost(elem(2,:)+1) - X_ghost(elem(1,:)+1);
t_y = Y_ghost(elem(2,:)+1) - Y_ghost(elem(1,:)+1);
ds = sqrt(t_x.^2 + t_y.^2);
% quiver(X_mid, Y_mid, t_x./ds, t_y./ds)
n_x = -t_y./ds;
n_y = t_x./ds;
quiver(X_mid, Y_mid, n_x, n_y)

end

