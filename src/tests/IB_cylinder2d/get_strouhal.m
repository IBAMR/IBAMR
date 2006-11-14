function [T,f,St] = get_strouhal(C_L)

  u_oo = 1.0;  % inflow velocity
  r = 0.15;    % radius of cylinder

  % compute the shedding times using linear interpolation
  T = [];
  for i = 1:size(C_L,1)-1
    if ((C_L(i,2) <= 0.0 && C_L(i+1,2) > 0.0) || (C_L(i,2) >= 0.0 && C_L(i+12) < 0.0))
      t0 = C_L(i  ,1);
      t1 = C_L(i+1,1);
      C0 = C_L(i  ,2);
      C1 = C_L(i+1,2);
      t_zero = (C0*t1-C1*t0)/(C0-C1);
      T = [T;t_zero];
    end
  end

  % compute the shedding frequency
  f = T(2:length(T)) - T(1:length(T)-1);

  % compute the Strouhal number
  St = 2.0./(u_oo*f/r);

  return;
