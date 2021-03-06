%% ---------------------------------------------------------------------
%%
%% Copyright (c) 2019 - 2019 by the IBAMR developers
%% All rights reserved.
%%
%% This file is part of IBAMR.
%%
%% IBAMR is free software and is distributed under the 3-clause BSD
%% license. The full text of the license can be found in the file
%% COPYRIGHT at the top level directory of IBAMR.
%%
%% ---------------------------------------------------------------------

i = 1;
figure(1); clf;
fileID = fopen(strcat('Sxx_component.curve'),'r');
while(true)
    tline = fgetl(fileID);
    if ~ischar(tline)
        break
    end
    A = str2num(tline);
    t = A(1);
    x = A(2:3:end);
    y = A(3:3:end);
    sxx = A(4:3:end);
    % find where y is non-zero
    idxs = find(y);
    x_new = x;
    x_new(idxs(1):idxs(end)) = sqrt(x(idxs(1):idxs(end)).^2+y(idxs(1):idxs(end)).^2);
    plot(x_new,sxx);
    title(['t = ', num2str(t)]);
    pause(0.01);
end
fclose(fileID);
