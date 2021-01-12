%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2014 - 2019 by the IBAMR developers
%% All rights reserved.
%%
%% This file is part of IBAMR.
%%
%% IBAMR is free software and is distributed under the 3-clause BSD
%% license. The full text of the license can be found in the file
%% COPYRIGHT at the top level directory of IBAMR.
%%
%% ---------------------------------------------------------------------

function n = L2normx(u)
   N1 = size(u,1);
   N2 = size(u,2);
   h = 1/N2;
   n = sqrt(0.5*sum(u(1,:).^2)*h^2 + sum(sum(u(2:N1-1,:).^2))*h^2 + 0.5*sum(u(N1,:).^2)*h^2);
