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

function n = L1normx(u)
   N1 = size(u,1);
   N2 = size(u,2);
   h = 1/N2;
   n = 0.5*sum(abs(u(1,:)))*h^2 + sum(sum(abs(u(2:N1-1,:))))*h^2 + 0.5*sum(abs(u(N1,:)))*h^2;
