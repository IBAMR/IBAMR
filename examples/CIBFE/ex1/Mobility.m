addpath('/home/donev/HPC/FluctHydro/IBAMR/ibamr_git/examples/CIBFE/ex1');
addpath('/home/fbalboa/sfw/linux/petsc/3.4.5/bin/matlab/');
MM = PetscBinaryRead('mobility_mat.dat');
MM=(MM+MM')/2; % Make it exactly symmetric

% Spectrum:
figure(1); clf;
semilogy(sort(eig(MM)), 'ko');

% Diagonal values -- effective hydrodynamic radius:
figure(2); 
plot(1./(6*pi*diagMM(1:3:end)), 'ro'); hold on;
plot(1./(6*pi*diagMM(2:3:end)), 'gs'); hold on;
plot(1./(6*pi*diagMM(3:3:end)), 'bd'); hold on;
