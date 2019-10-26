function [X, ds]  = generate_leaf2d(N)

  L = 16.0;
  dx = L/N
  ds = 1*dx;

  X0 = 0;
  Y0 = 0;
  l_x = 1;
  l_y = 1;
  n_x = ceil(l_x/ds);
  if (mod(n_x,2) == 0)
    n_x = n_x + 1;
  end
  n_y = ceil(l_y/ds);
  if (mod(n_y,2) == 0)
    n_y = n_y + 1;
  end
  ds_x = l_x/(n_x-1);
  ds_y = l_y/(n_y-1);

vertex_fid = fopen(['stalk2d_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', n_x);
for n = 1:n_x
     fprintf(vertex_fid, '%1.16e %1.16e\n', X0 + (n-1)*ds_x, Y0);
end 
fclose(vertex_fid);

kappa = 0.0;
spring_fid = fopen(['stalk2d_' num2str(N) '.spring'], 'w');
fprintf(spring_fid, '%d\n', n_x-1);
for n = 1:n_x-1
     fprintf(spring_fid, '%6d %6d %1.16e %1.16e\n', n-1, n, kappa, ds_x);
end
fclose(spring_fid);

kappa = 0.0;
beam_fid = fopen(['stalk2d_' num2str(N) '.beam'], 'w');
fprintf(beam_fid, '%d\n', n_x-2);
for n = 1:n_x-2
     fprintf(beam_fid, '%6d %6d %6d %1.16e\n', n-1, n, n+1, kappa);
end
fclose(beam_fid);

target_fid = fopen(['stalk2d_' num2str(N) '.target'], 'w');
fprintf(target_fid, '%d\n', 1);
fprintf(target_fid, '%6d %1.16e\n', 0, kappa);
fclose(target_fid);

vertex_fid = fopen(['leaf2d_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', n_y);
for n = 1:n_y
     fprintf(vertex_fid, '%1.16e %1.16e\n', X0 + l_x, Y0 - 0.5*l_y + (n-1)*ds_y);
end
fclose(vertex_fid);

kappa = 0.0;
spring_fid = fopen(['leaf2d_' num2str(N) '.spring'], 'w');
fprintf(spring_fid, '%d\n', n_y-1);
for n = 1:n_y-1
     fprintf(spring_fid, '%6d %6d %1.16e %1.16e\n', n-1, n, kappa, ds_y);
end
fclose(spring_fid);

kappa = 0.0;
beam_fid = fopen(['leaf2d_' num2str(N) '.beam'], 'w');
fprintf(beam_fid, '%d\n', n_y-2);
for n = 1:n_y-2
     fprintf(beam_fid, '%6d %6d %6d %1.16e\n', n-1, n, n+1, kappa);
end
fclose(beam_fid);

