%
% initialize the IB points for a thick ring of points and write 
% down the corresponding vertex and spring files
%
function [X,ds, ds2]  = generate_shell2d(N)
%% geometric parameters for the circle
  %
  xc = 0.5;
% center of the circle yc = 0.5;
% r = 0.25;
% radius of the circle w = 0.0625;   % width of the band
  
  % define the number of points in each direction
  %  these formula are for r=1/4, w=1/16 and give spacing
  %  of ds~2dx/3
  %
  M1 = 19/8*N;
M2 = 3 / 32 * N + 1;

% define ds based on the circumfrential direction % ds = 2 * pi* r / M1 ds2 = w / (M2 - 1);
theta = (linspace(0, 2 * pi, M1 + 1))'; theta = theta(1 : end - 1);

% build the data strcutres X = [];
  for
      k = 1 : M2

                  rk = r + ds2 * (k - 1);
  Xk = [ xc + rk * cos(theta), yc + rk* sin(theta) ];
  X = [X; Xk];
  end

      num_nodes = size(X, 1);
  num_nodes_ring = num_nodes / M2;
  num_springs = size(X, 1);

  % % Step 1 : Write out the vertex information vertex_fid = fopen(['shell2d_' num2str(N)'.vertex'], 'w');
  fprintf(vertex_fid, '%d\n', num_nodes);

for
    node = 1 : num_nodes fprintf(vertex_fid, '%1.16e %1.16e\n', X(node, 1), X(node, 2));
end

    % % Step 2 : Write out the spring information spring_fid = fopen(['shell2d_' num2str(N)'.spring'], 'w');
fprintf(spring_fid, '%d\n', num_springs);

rest_length = 0.0;
Kappa = 0.0;
for
    ring = 1 : M2 lag_begin = (ring - 1) * num_nodes_ring;
lag_end = lag_begin + (num_nodes_ring - 1);
    for
        node = lag_begin : lag_end -
                           1 fprintf(spring_fid, '%6d %6d %1.16e %1.16e\n', node, node + 1, ... Kappa, rest_length);
    end fprintf(spring_fid, '%6d %6d %1.16e %1.16e\n', lag_end, lag_begin, ... Kappa, rest_length);
    end
