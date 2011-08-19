%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem parameters
L = 10.0;         % length of computational domain (cm)
N = 64;           % number of Cartesian grid points
h = L/N;          % Cartesian grid spacing
rho = 1.0;        % fluid density (g/cm^3)
mu = 0.01;        % fluid viscosity (g/cm*s)

l0 = 5.0;         % length of rod (cm)
nr = 65;          % number of Lagrangian mesh nodes
ds = l0/(nr-1);   % Lagrangian meshwidth

dt = 0.01;        % timestep size (s)

a = 0.3;
a1 = a;
a2 = a;
a3 = (2/3)*a;

b = 54;
b1 = b;
b2 = b;
b3 = b;

kappa1 = 0.0;
kappa2 = 0.0;
tau = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertex_fid = fopen(['curve3d_' num2str(N) '.vertex'], 'w');
rod_fid = fopen(['curve3d_' num2str(N) '.rod'], 'w');
director_fid = fopen(['curve3d_' num2str(N) '.director'], 'w');

fprintf(vertex_fid, '%d\n', nr);
fprintf(rod_fid, '%d\n', nr-1);
fprintf(director_fid, '%d\n', nr);

for r = 0:nr-1
  X = [0 0 l0*(r/(nr-1) - 0.5)];

  fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', X(1), X(2), X(3));

  r_curr = r;
  r_next = r+1;

  if (r_next < nr)
    fprintf(rod_fid, ['%6d %6d %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n'], ...
            r_curr, r_next, ds, a1, a2, a3, b1, b2, b3, kappa1, kappa2, tau);
  end

  D1 = [1 0 0];
  D2 = [0 1 0];
  D3 = [0 0 1];

  fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D1(1), D1(2), D1(3));
  fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D2(1), D2(2), D2(3));
  fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D3(1), D3(2), D3(3));
end %for

fclose(vertex_fid);
fclose(rod_fid);
fclose(director_fid);

% Anchor the first point in space.
anchor_fid = fopen(['curve3d_' num2str(N) '.anchor'], 'w');
fprintf(anchor_fid, '%d\n', 1);
fprintf(anchor_fid, '%6d\n', 0);
fclose(anchor_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
