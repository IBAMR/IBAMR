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

 A = [1,0,0,0;
     1,1/3, 1/9, 1/27;
	 1, 2/3, 4/9, 8/27;
	 1, 1  ,  1 ,  1];
 
 K = [1.29;0.52;5.43;4.28];
 
 coefs = inv(A)*K
  