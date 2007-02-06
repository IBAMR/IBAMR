function St = compute_St(C_L,dt)
  r    = 0.15;  % cylinder radius
  D    = 2.0*r; % cylinder diameter
  u_oo = 1.00;  % inflow velocity

  n = length(C_L);
  Fs = 1/dt;
  f = (0:floor(n/2))*(Fs/n);

  p = abs(fft(C_L));
  p = p(1:floor(n/2)+1);

  [max_p , i] = max(p);
  St = f(i)*D/u_oo;

  return;
