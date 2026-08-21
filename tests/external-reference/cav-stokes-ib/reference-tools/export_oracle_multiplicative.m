function export_oracle_multiplicative()
% Export one pinned-sandbox multiplicative-CAV level and native control.

oracle_root = getenv('CAV_ORACLE_ROOT');
output_dir = getenv('CAV_ORACLE_OUTPUT');
physical_K = str2double(getenv('CAV_PHYSICAL_K'));
if isempty(oracle_root) || isempty(output_dir) || ~isfinite(physical_K) || physical_K <= 0
    error(['CAV_ORACLE_ROOT, CAV_ORACLE_OUTPUT, and a positive ' ...
           'CAV_PHYSICAL_K are required.']);
end

restoredefaultpath();
addpath(oracle_root);
setup_project_paths(oracle_root);

N = 8;
n = N^2;
h = 1/N;
dt = 0.005;
rho = 1;
mu = 0.01;
ds = h;
center = [0.5, 0.5];
radius = [0.2, 0.2];
marker_count = max(3, floor(2*pi*radius(1)/ds));
theta = (2*pi)*(0:(marker_count-1))'/marker_count;
X = [center(1)+radius(1)*cos(theta), center(2)+radius(2)*sin(theta)];

[A, ib_data] = build_saddle_point_matrix_sparse( ...
    X, N, h, h, ds, physical_K, dt, rho, mu, 'spring');
if isfield(ib_data, 'IB_uu')
    ib_velocity = [ib_data.IB_uu, ib_data.IB_uv; ...
                   ib_data.IB_vu, ib_data.IB_vv];
else
    Z = sparse(n, n);
    ib_velocity = [ib_data.IB_u, Z; Z, ib_data.IB_v];
end
E_h = blkdiag(-dt*ib_velocity, sparse(n, n));
[smoother_data, ~] = precompute_smoother_data( ...
    A, N, 1, [], 'targeted_ib', ib_velocity);

% This state is intentionally deterministic and independent of the candidate
% export. Its assembled RHS gives both patch sequences a nontrivial control.
probe = zeros(3*n, 1);
for i = 0:N-1
    for j = 0:N-1
        k = j+N*i+1;
        probe(k) = sin(2*pi*(j+0.5)/N)+0.25*cos(4*pi*i/N);
        probe(n+k) = -sin(2*pi*(i+0.5)/N)+0.125*cos(2*pi*j/N);
        probe(2*n+k) = cos(2*pi*i/N)*sin(2*pi*j/N);
    end
end
probe(2*n+1:end) = probe(2*n+1:end)-mean(probe(2*n+1:end));
rhs = A*probe;
[correction, incremental_residual] = coupling_based_smoother_precomp( ...
    A, rhs, zeros(3*n, 1), N, 1, smoother_data, []);
fresh_residual = rhs-A*correction;

patch0 = smoother_data.coupled_dofs{1};
patch0_matrix = A(patch0, patch0);
patch0_rhs = rhs(patch0);
patch0_correction = patch0_matrix\patch0_rhs;

write_matrix_market(fullfile(output_dir, 'oracle_A.mtx'), A);
write_matrix_market(fullfile(output_dir, 'oracle_E_h.mtx'), E_h);
write_vector_market(fullfile(output_dir, 'oracle_probe.mtx'), probe);
write_vector_market(fullfile(output_dir, 'oracle_rhs.mtx'), rhs);
write_vector_market(fullfile(output_dir, 'oracle_smoother_correction.mtx'), correction);
write_vector_market(fullfile(output_dir, 'oracle_smoother_incremental_residual.mtx'), incremental_residual);
write_vector_market(fullfile(output_dir, 'oracle_smoother_fresh_residual.mtx'), fresh_residual);
write_matrix_market(fullfile(output_dir, 'oracle_patch0_A.mtx'), patch0_matrix);
write_vector_market(fullfile(output_dir, 'oracle_patch0_rhs.mtx'), patch0_rhs);
write_vector_market(fullfile(output_dir, 'oracle_patch0_correction.mtx'), patch0_correction);
write_dof_map(fullfile(output_dir, 'oracle_dof_map.txt'), N);
write_index_list(fullfile(output_dir, 'oracle_pressure_seeds.txt'), ...
                 2*n+smoother_data.seed_dofs-1);
write_index_sets(fullfile(output_dir, 'oracle_patches.txt'), ...
                 smoother_data.coupled_dofs);

metadata = struct();
metadata.schema = 'sandbox-cav-multiplicative-oracle-v2';
metadata.oracle_sha = '5b77344db6746269f8c77695c99e9043907ba74b';
metadata.dimension = int32(2);
metadata.base_grid_cells_per_axis = int32(4);
metadata.refinement_ratio = int32(2);
metadata.finest_grid_cells_per_axis = int32(N);
metadata.dt = dt;
metadata.rho = rho;
metadata.mu = mu;
metadata.lagrangian_spacing = ds;
metadata.physical_stiffness_K = physical_K;
metadata.marker_count = int32(marker_count);
metadata.structure_center = center;
metadata.structure_radius = radius;
metadata.delta_function = 'IB_4';
metadata.patch_strategy = 'targeted_ib';
metadata.seed_stride = int32(1);
metadata.traversal_order = 'pressure linear index, y-fastest';
metadata.composition = 'multiplicative';
metadata.pressure_equation = 'div';
metadata.pressure_gauge = 'zero-mean correction';
metadata.patch_count = int32(numel(smoother_data.coupled_dofs));
metadata.standard_patch_count = int32(smoother_data.standard_patch_count);
metadata.enlarged_patch_count = int32(smoother_data.enlarged_patch_count);
write_json(fullfile(output_dir, 'oracle_metadata.json'), metadata);

fprintf('oracle_physical_stiffness_K=%.17g\n', physical_K);
fprintf('oracle_patch_count=%d\n', numel(smoother_data.coupled_dofs));
fprintf('oracle_enlarged_patch_count=%d\n', smoother_data.enlarged_patch_count);
end

function write_dof_map(filename, N)
n = N^2;
fid = open_output(filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'ibamr-cav-global-dof-map-v1 2 %d\n', 3*n);
for axis = 0:1
    for i = 0:N-1
        for j = 0:N-1
            block_dof = j+N*i;
            dof = axis*n+block_dof;
            if axis == 0
                x = i/N;
                y = (j+0.5)/N;
            else
                x = (i+0.5)/N;
                y = j/N;
            end
            fprintf(fid, '%d V %d %d %d %.17e %.17e\n', ...
                    dof, axis, i, j, x, y);
        end
    end
end
for i = 0:N-1
    for j = 0:N-1
        dof = 2*n+j+N*i;
        fprintf(fid, '%d P -1 %d %d %.17e %.17e\n', ...
                dof, i, j, (i+0.5)/N, (j+0.5)/N);
    end
end
end

function write_index_list(filename, values)
fid = open_output(filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'ibamr-cav-index-list-v1 %d\n', numel(values));
fprintf(fid, '%d\n', values(:));
end

function write_index_sets(filename, sets)
fid = open_output(filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'ibamr-cav-index-sets-v1 %d\n', numel(sets));
for ordinal = 1:numel(sets)
    values = sort(sets{ordinal}(:)-1);
    fprintf(fid, '%d %d', ordinal-1, numel(values));
    fprintf(fid, ' %d', values);
    fprintf(fid, '\n');
end
end

function write_matrix_market(filename, A)
[rows, columns, values] = find(A);
fid = open_output(filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%%%%MatrixMarket matrix coordinate real general\n');
fprintf(fid, '%d %d %d\n', size(A, 1), size(A, 2), numel(values));
for k = 1:numel(values)
    fprintf(fid, '%d %d %.17e\n', rows(k), columns(k), values(k));
end
end

function write_vector_market(filename, vector)
fid = open_output(filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%%%%MatrixMarket matrix array real general\n');
fprintf(fid, '%d 1\n', numel(vector));
fprintf(fid, '%.17e\n', vector(:));
end

function write_json(filename, value)
fid = open_output(filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(value, PrettyPrint=true));
end

function fid = open_output(filename)
fid = fopen(filename, 'w');
if fid < 0
    error('Unable to open %s for writing.', filename);
end
end
