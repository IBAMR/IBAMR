function f = f_i(l)
  mu    = 0.034;  % MPa
  k1    = 4.34;   % MPa
  k2    = 13.32;  % dimensionless
  kappa = 0.20;   % dimensionless

  FF = eye(3,3);
  FF_0 = l*FF;
  CC_0 = FF_0'*FF_0;
  I1_0 = trace(CC_0);
  I4_star = l^2;
  f = k1*(I4_star-1.0)*exp(k2*(I4_star-1.0)^2)*l^2;