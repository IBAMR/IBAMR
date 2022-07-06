%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2022 - 2022 by the IBAMR developers
%% All rights reserved.
%%
%% This file is part of IBAMR.
%%
%% IBAMR is free software and is distributed under the 3-clause BSD
%% license. The full text of the license can be found in the file
%% COPYRIGHT at the top level directory of IBAMR.
%%
%% ---------------------------------------------------------------------

% Problem parameters
L = 0.08;                 % computational domain linear dimension (cm)
ASPECT_RATIO_X = 1.0;     % ratio of computational domain length and characteristic linear dimension (dimensionless)
ASPECT_RATIO_Y = 1.0;     % ratio of computational domain width  and characteristic linear dimension (dimensionless)
ASPECT_RATIO_Z = 0.025;    % ratio of computational domain height and characteristic linear dimension (dimensionless)

post_density = 7.0e4;     % posts per unit area (cm^{-2}) (updated Jun 30 2021)
post_length = 6.0e-3;
post_deflection_radius = 2.95e-3; %(cm)
slanted_post_height = sqrt(post_length^2 - post_deflection_radius^2);% post height (cm) (updated Jun 30 2021)


n_posts = 1; %Just one post for visualization purposes

NFINEST = 128;            % number of grid cells on finest grid level
h = L/NFINEST;            % Cartesian grid spacing

n_ib_post = 30;  % number of IB points per post
n_ib = n_ib_post * n_posts;                % total number of IB points
dX = post_length / (n_ib_post-1);          % IB point spacing along the post

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertex_fid = fopen(['posts_' num2str(NFINEST) '.vertex'], 'w');


fprintf(vertex_fid, '%d\n', n_posts * n_ib_post);


for p = 0:n_posts-1

  % vertices:
  X(1) = 0.08/2;
X(2) = 0.08/2; %(x,y) coordinates which place the post in the center of the domain
  for l = 0:n_ib_post-1
     fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', (post_deflection_radius/post_length)*(l*dX) + X(1), X(2), (slanted_post_height/post_length)*(l*dX));
  end

end

fclose(vertex_fid);
