%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters
L = 10.0;         % length of computational domain (cm)
N = 64;           % number of Cartesian grid points
h = L/N;          % Cartesian grid spacing
rho = 1.0;        % fluid density (g/cm^3)
mu = 0.01;        % fluid viscosity (g/cm*s)

r0 = 2.5;         % radius of unstressed rod (cm)
nr = 200;         % number of Lagrangian mesh nodes
ds = 2*pi*r0/nr;  % Lagrangian meshwidth

dt = 0.01;        % timestep size (s)

p = 5;            % number of twists
epsilon = 0.001;  % perturbation parameter

a = 0.3;
a3 = (2/3)*a;
b = 54;
beta = asin(-(a3*p)/(b*r0^2+a3-a));
r1 = r0*cos(beta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertex_fid = fopen(['curve3d_' num2str(N) '.vertex'], 'w');
beam_fid = fopen(['curve3d_' num2str(N) '.beam'], 'w');
director_fid = fopen(['curve3d_' num2str(N) '.director'], 'w');

fprintf(vertex_fid, '%d\n', nr);
fprintf(beam_fid, '%d\n', nr);
fprintf(director_fid, '%d\n', nr);

for r = 0:nr-1
  theta = 2*pi*r/nr;

  R = [cos(theta) , sin(theta) , 0];
  Theta = [-sin(theta), cos(theta) , 0];
  Z = [0 , 0 , 1];

  X = r1*R;

  fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', X(1), X(2), X(3));

  r_next = mod(r+1,nr);
  r_prev = mod(r-1,nr);

  fprintf(beam_fid, '%6d %6d %6d %1.16e\n', r_prev, r, r_next, 0.0);

  D3 =  cos(beta)*Theta + sin(beta)*Z;
  E  = -sin(beta)*Theta + cos(beta)*Z;
  D1 =  cos(p*theta + epsilon*sin(theta))*E + sin(p*theta + epsilon*sin(theta))*R;
  D2 = -sin(p*theta + epsilon*sin(theta))*E + cos(p*theta + epsilon*sin(theta))*R;

  fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D1(1), D1(2), D1(3));
  fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D2(1), D2(2), D2(3));
  fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D3(1), D3(2), D3(3));
end %for

fclose(vertex_fid);
fclose(beam_fid);
fclose(director_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
