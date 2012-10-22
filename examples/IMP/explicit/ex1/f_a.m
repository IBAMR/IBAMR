function f = f_a(l)
  mu    = 0.020;  % MPa
  k1    = 0.39;   % MPa
  k2    = 6.79;   % dimensionless
  kappa = 0.23;   % dimensionless

  FF = eye(3,3);
  FF_0 = l*FF;
  CC_0 = FF_0'*FF_0;
  I1_0 = trace(CC_0);
  I4_star = l^2;
  f = k1*(I4_star-1.0)*exp(k2*(I4_star-1.0)^2)*l^2;