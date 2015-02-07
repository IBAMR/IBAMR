%clear all
d=2;
addpath('/home/donev/HPC/FluctHydro/IBAMR/ibamr_git/examples/CIBFE/ex1');
addpath('/home/fbalboa/sfw/linux/petsc/3.4.5/bin/matlab/');
MM = PetscBinaryRead('mobility_mat.dat');
MM=(MM+MM')/2; % Make it exactly symmetric
diagMM=diag(MM);
ratio = max(diagMM)/min(diagMM)
eigMM=eig(MM);
cond_num = max(eigMM)/min(eigMM)

% Spectrum:
figure(1); clf;
semilogy(sort(eigMM/max(eigMM)), 'ko');

% Diagonal values -- effective hydrodynamic radius:
figure(2); clf;
plot(1./(6*pi*diagMM(1:d:end)), 'ro'); hold on;
plot(1./(6*pi*diagMM(2:d:end)), 'gs'); hold on;

% Now also read and include mass matrix:
% --------------------
mass = PetscBinaryRead('mass_mat.dat');
invmass = inv(mass); % Form explicit inverse here for simplicity
%figure(3); clf; semilogy(sort(eig(mass)), 'ko');

scaledMM = invmass * MM * invmass;

MM=(scaledMM+scaledMM')/2; % Make it exactly symmetric
diagMM=diag(MM);
ratio = max(diagMM)/min(diagMM)
eigMM=eig(MM);
cond_num = max(eigMM)/min(eigMM)

% Spectrum:
figure(1); hold on;
semilogy(sort(eigMM/max(eigMM)), 'rs');
legend('No mass','With mass');

% Diagonal values -- effective hydrodynamic radius:
figure(3); clf;
plot(1./(6*pi*diagMM(1:d:end)), 'ro'); hold on;
plot(1./(6*pi*diagMM(2:d:end)), 'gs'); hold on;

