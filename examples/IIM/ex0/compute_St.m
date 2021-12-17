%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2014 - 2014 by the IBAMR developers
%% All rights reserved.
%%
%% This file is part of IBAMR.
%%
%% IBAMR is free software and is distributed under the 3-clause BSD
%% license. The full text of the license can be found in the file
%% COPYRIGHT at the top level directory of IBAMR.
%%
%% ---------------------------------------------------------------------

function [St_fft St_interp] = compute_St(C_L,dt)
  r    = 0.5;   % cylinder radius
  D    = 2.0*r; % cylinder diameter
  u_oo = 1.00;  % inflow velocity

  % determine the Strouhal number via discrete Fourier analysis
  n = length(C_L);
  Fs = 1/dt;
  f = (0:floor(n/2))*(Fs/n);

  p = abs(fft(C_L));
  p = p(1:floor(n/2)+1);

  [max_p , i] = max(p);
  St_fft = f(i)*D/u_oo;

  plot(f,p,'-',f(i),max_p,'o');
  xlabel('frequency');
  ylabel('power');
  title('discrete Fourier analysis of C_L');

  % determine the Strouhal number by linear interpolation
  a = 1;
  b = ones(4,1)/4;

  C_L = filter(b,a,C_L);

  t_previous = 0.0;
  t_current = 0.0;
  St_interp = [];
  for i = 1:length(C_L)-1
    if (C_L(i+1) >= 0.0 & C_L(i) < 0.0)
      t_previous = t_current;

      t_start = (i-1)*dt;   % time at the beginning of the interval
      t_end = i*dt;         % time at the end       of the interval

      C_L_start = C_L(i-1); % C_L at the beginning of the interval
      C_L_end = C_L(i);     % C_L at the end       of the interval

      % the time when the interpolated C_L crosses zero
      t_current = (C_L_start*t_end - C_L_end*t_start)/(C_L_start - C_L_end);

      St_interp = [St_interp (1.0/(t_current-t_previous))*D/u_oo];
    end %if
  end %for

  St_interp = St_interp(length(St_interp)-20:length(St_interp));
  if (var(St_interp) > sqrt(eps))
    warning('high variance in the estimated St');
  end %if
  St_interp = mean(St_interp);

  return;
