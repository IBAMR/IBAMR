function replay_n0_common_arithmetic()
% Compare N0 construction and replay both patch sequences in one runtime.

candidate_dir = getenv('CAV_N0_CANDIDATE_OUTPUT');
oracle_dir = getenv('CAV_N0_ORACLE_OUTPUT');
report_path = getenv('CAV_N0_REPORT');
if isempty(candidate_dir) || isempty(oracle_dir) || isempty(report_path)
    error(['CAV_N0_CANDIDATE_OUTPUT, CAV_N0_ORACLE_OUTPUT, and ' ...
           'CAV_N0_REPORT are required.']);
end

target = 1e-13;
candidate_map = read_dof_map(fullfile(candidate_dir, ...
    'cav_live_export_contract_dof_map.txt'));
oracle_map = read_dof_map(fullfile(oracle_dir, 'oracle_dof_map.txt'));
[mapping, mapping_report] = build_mapping(candidate_map, oracle_map);
pressure = strcmp(oracle_map.kind, 'P');

candidate_A_raw = read_matrix_market(fullfile(candidate_dir, ...
    'cav_live_export_contract_A.mtx'));
candidate_E_raw = read_matrix_market(fullfile(candidate_dir, ...
    'cav_live_export_contract_E_h.mtx'));
oracle_A = read_matrix_market(fullfile(oracle_dir, 'oracle_A.mtx'));
oracle_E = read_matrix_market(fullfile(oracle_dir, 'oracle_E_h.mtx'));
candidate_A = map_matrix(candidate_A_raw, mapping, pressure, true);
candidate_E = map_matrix(candidate_E_raw, mapping, pressure, false);

operator_comparison = compare_matrices(candidate_A, oracle_A);
elasticity_comparison = compare_matrices(candidate_E, oracle_E);

candidate_seeds_raw = read_index_list(fullfile(candidate_dir, ...
    'cav_live_export_contract_pressure_seeds.txt'));
candidate_patches_raw = read_index_sets(fullfile(candidate_dir, ...
    'cav_live_export_contract_patches.txt'));
oracle_seeds = read_index_list(fullfile(oracle_dir, 'oracle_pressure_seeds.txt'))+1;
oracle_patches = one_based_sets(read_index_sets(fullfile(oracle_dir, ...
    'oracle_patches.txt')));
candidate_seeds = mapping(candidate_seeds_raw+1);
candidate_patches = map_sets(candidate_patches_raw, mapping);
patch_comparison = compare_patch_data( ...
    candidate_seeds, candidate_patches, oracle_seeds, oracle_patches);

candidate_rhs = map_vector(read_vector_market(fullfile(candidate_dir, ...
    'cav_fac_stage0_pre-smooth-input_level1_rhs.mtx')), mapping, pressure, true);
candidate_initial = map_vector(read_vector_market(fullfile(candidate_dir, ...
    'cav_fac_stage0_pre-smooth-input_level1_solution.mtx')), mapping, pressure, false);
candidate_actual = map_vector(read_vector_market(fullfile(candidate_dir, ...
    'cav_fac_stage1_pre-smooth-output_level1_solution.mtx')), mapping, pressure, false);
oracle_rhs = read_vector_market(fullfile(oracle_dir, 'oracle_rhs.mtx'));

[live_candidate_w, live_candidate_r] = apply_sweep( ...
    candidate_A, candidate_rhs, candidate_initial, candidate_patches, pressure);
[live_oracle_w, live_oracle_r] = apply_sweep( ...
    candidate_A, candidate_rhs, candidate_initial, oracle_patches, pressure);
[oracle_candidate_w, oracle_candidate_r] = apply_sweep( ...
    oracle_A, oracle_rhs, zeros(size(oracle_rhs)), candidate_patches, pressure);
[oracle_oracle_w, oracle_oracle_r, oracle_patch0] = apply_sweep( ...
    oracle_A, oracle_rhs, zeros(size(oracle_rhs)), oracle_patches, pressure);

live_common_replay = paired_result( ...
    live_candidate_w, live_candidate_r, live_oracle_w, live_oracle_r, pressure, target);
oracle_common_replay = paired_result( ...
    oracle_candidate_w, oracle_candidate_r, oracle_oracle_w, oracle_oracle_r, pressure, target);
actual_live_control = struct( ...
    'correction', compare_vectors_gauged(live_candidate_w, candidate_actual, pressure), ...
    'fresh_residual', compare_vectors(candidate_rhs-candidate_A*live_candidate_w, ...
                                      candidate_rhs-candidate_A*candidate_actual));
actual_live_control.within_fixed_target = ...
    actual_live_control.correction.e_inf <= target && ...
    actual_live_control.fresh_residual.e_inf <= target;
actual_live_control.classification = ...
    'reported cross-runtime backend/arithmetic diagnostic; not a common-arithmetic gate';

oracle_actual = read_vector_market(fullfile(oracle_dir, ...
    'oracle_smoother_correction.mtx'));
oracle_actual_r = read_vector_market(fullfile(oracle_dir, ...
    'oracle_smoother_fresh_residual.mtx'));
oracle_pinned_control = struct( ...
    'correction', compare_vectors_gauged(oracle_oracle_w, oracle_actual, pressure), ...
    'fresh_residual', compare_vectors(oracle_oracle_r, oracle_actual_r));
oracle_pinned_control.within_fixed_target = ...
    oracle_pinned_control.correction.e_inf <= target && ...
    oracle_pinned_control.fresh_residual.e_inf <= target;
oracle_pinned_control.classification = ...
    'reported cross-runtime backend/arithmetic diagnostic; not a common-arithmetic gate';

candidate_local_control = check_candidate_local_trace( ...
    candidate_dir, candidate_A_raw, candidate_patches_raw{1});
oracle_local_control = struct( ...
    'matrix', compare_matrices(oracle_patch0.A, read_matrix_market( ...
        fullfile(oracle_dir, 'oracle_patch0_A.mtx'))), ...
    'rhs', compare_vectors(oracle_patch0.rhs, read_vector_market( ...
        fullfile(oracle_dir, 'oracle_patch0_rhs.mtx'))), ...
    'correction', compare_vectors(oracle_patch0.delta, read_vector_market( ...
        fullfile(oracle_dir, 'oracle_patch0_correction.mtx'))));
oracle_local_control.pass = oracle_local_control.matrix.e_inf <= target && ...
                            oracle_local_control.rhs.e_inf <= target && ...
                            oracle_local_control.correction.e_inf <= target;

mutations = run_mutations(candidate_A_raw, candidate_A, oracle_A, candidate_map, ...
                          oracle_map, mapping, pressure, candidate_patches, ...
                          oracle_patches, live_candidate_w, target);

report = struct();
report.schema = 'cav-n0-paired-common-arithmetic-v1';
report.scope = ['supported 4x4-to-8x8 two-level 2D hierarchy; finest-level ' ...
                'RELAXED pressure-seed multiplicative CAV'];
report.physical_stiffness_K = 1;
report.fixed_target = target;
report.error_formula = ['E_inf=max(abs(x-y))/max(1,max(abs(x)),' ...
                        'max(abs(y))); matrices use the same formula over all entries'];
report.runtime = struct('name', 'MATLAB', 'version', version, 'computer', computer);
report.mapping = mapping_report;
report.operator = operator_comparison;
report.elasticity = elasticity_comparison;
report.patches = patch_comparison;
report.live_operator_common_arithmetic = live_common_replay;
report.oracle_operator_common_arithmetic = oracle_common_replay;
report.actual_live_smoother_diagnostic_not_gate = actual_live_control;
report.pinned_oracle_smoother_diagnostic_not_gate = oracle_pinned_control;
report.candidate_selected_local_solve = candidate_local_control;
report.oracle_selected_local_solve = oracle_local_control;
report.mutations = mutations;
report.pass = mapping_report.pass && operator_comparison.structure_exact && ...
    elasticity_comparison.structure_exact && patch_comparison.pass && ...
    live_common_replay.pass && oracle_common_replay.pass && ...
    candidate_local_control.pass && oracle_local_control.pass && mutations.pass;

write_json(report_path, report);
fprintf('n0_operator_e_inf=%.17e\n', operator_comparison.e_inf);
fprintf('n0_elasticity_e_inf=%.17e\n', elasticity_comparison.e_inf);
fprintf('n0_patch_structure_exact=%d\n', patch_comparison.pass);
fprintf('n0_paired_replay_pass=%d\n', report.pass);
if ~report.pass
    error('N0 paired common-arithmetic replay failed; see %s.', report_path);
end
end

function [mapping, result] = build_mapping(source, target)
if source.count ~= target.count || source.dimension ~= target.dimension
    error('Candidate and oracle DOF-map dimensions differ.');
end
target_by_key = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:target.count
    key = dof_key(target, k);
    if isKey(target_by_key, key)
        error('Duplicate oracle DOF key %s.', key);
    end
    target_by_key(key) = k;
end
mapping = zeros(source.count, 1);
max_position_error = 0;
for k = 1:source.count
    key = dof_key(source, k);
    if ~isKey(target_by_key, key)
        error('Candidate DOF key is absent from oracle map: %s.', key);
    end
    mapping(k) = target_by_key(key);
    max_position_error = max(max_position_error, ...
        max(abs(source.position(k, :)-target.position(mapping(k), :))));
end
bijective = numel(unique(mapping)) == source.count && ...
            isequal(sort(mapping), (1:source.count)');
result = struct('count_exact', source.count == target.count, ...
                'logical_keys_bijective', bijective, ...
                'max_position_error', max_position_error, ...
                'position_target', 1e-14, ...
                'pass', bijective && max_position_error <= 1e-14);
end

function key = dof_key(map, k)
key = sprintf('%s:%d:%d:%d', map.kind{k}, map.axis(k), map.i(k), map.j(k));
end

function mapped = map_matrix(raw, mapping, pressure, equation_sign)
mapped = sparse(numel(mapping), numel(mapping));
mapped(mapping, mapping) = raw;
if equation_sign
    mapped(pressure, :) = -mapped(pressure, :);
end
end

function mapped = map_vector(raw, mapping, pressure, equation_sign)
mapped = zeros(size(raw));
mapped(mapping) = raw;
if equation_sign
    mapped(pressure) = -mapped(pressure);
end
end

function mapped = map_sets(raw_sets, mapping)
mapped = cell(size(raw_sets));
for k = 1:numel(raw_sets)
    mapped{k} = sort(mapping(raw_sets{k}+1));
end
end

function sets = one_based_sets(raw_sets)
sets = cell(size(raw_sets));
for k = 1:numel(raw_sets)
    sets{k} = sort(raw_sets{k}+1);
end
end

function result = compare_patch_data(left_seeds, left, right_seeds, right)
result = struct('patch_count_exact', numel(left) == numel(right), ...
                'seed_order_exact', isequal(left_seeds(:), right_seeds(:)), ...
                'membership_and_order_exact', false, ...
                'first_mismatch_ordinal', int32(-1));
if result.patch_count_exact
    result.membership_and_order_exact = true;
    for k = 1:numel(left)
        if ~isequal(left{k}(:), right{k}(:))
            result.membership_and_order_exact = false;
            result.first_mismatch_ordinal = int32(k-1);
            break;
        end
    end
end
result.pass = result.patch_count_exact && result.seed_order_exact && ...
              result.membership_and_order_exact;
end

function [w, fresh_residual, first_patch] = apply_sweep(A, b, w, patches, pressure)
r = b-A*w;
first_patch = struct();
for patch = 1:numel(patches)
    dofs = patches{patch};
    local_A = full(A(dofs, dofs));
    local_rhs = r(dofs);
    delta = local_A\local_rhs;
    if patch == 1
        first_patch = struct('A', local_A, 'rhs', local_rhs, 'delta', delta);
    end
    w(dofs) = w(dofs)+delta;
    r = r-A(:, dofs)*delta;
end
w(pressure) = w(pressure)-mean(w(pressure));
fresh_residual = b-A*w;
end

function result = paired_result(left_w, left_r, right_w, right_r, pressure, target)
result = struct('correction', compare_vectors_gauged(left_w, right_w, pressure), ...
                'fresh_residual', compare_vectors(left_r, right_r));
result.pass = result.correction.e_inf <= target && result.fresh_residual.e_inf <= target;
end

function result = check_candidate_local_trace(candidate_dir, A, patch0_zero_based)
prefix = fullfile(candidate_dir, 'cav_fgmres_local_sweep0_patch0');
local_A = read_matrix_market([prefix '_A.mtx']);
local_rhs = read_vector_market([prefix '_rhs.mtx']);
local_correction = read_vector_market([prefix '_correction.mtx']);
pre_residual = read_vector_market([prefix '_pre_update_residual.mtx']);
dofs = patch0_zero_based+1;
extracted_A = A(dofs, dofs);
extracted_rhs = pre_residual(dofs);
backward_numerator = norm(local_A*local_correction-local_rhs, inf);
backward_denominator = max(1, norm(local_A, inf)*norm(local_correction, inf)+norm(local_rhs, inf));
result = struct( ...
    'matrix', compare_matrices(extracted_A, local_A), ...
    'restricted_rhs', compare_vectors(extracted_rhs, local_rhs), ...
    'direct_solve', compare_vectors(local_A\local_rhs, local_correction), ...
    'backward_error', backward_numerator/backward_denominator);
result.pass = result.matrix.e_inf <= 1e-13 && ...
              result.restricted_rhs.e_inf <= 1e-13 && ...
              result.direct_solve.e_inf <= 1e-13 && result.backward_error <= 1e-13;
end

function result = run_mutations(raw_A, mapped_A, oracle_A, candidate_map, oracle_map, ...
                                mapping, pressure, patches, oracle_patches, state, target)
count = numel(mapping);
permutation = [(2:count) 1]';
permuted_map = permute_dof_map(candidate_map, permutation);
[declared_mapping, declared_report] = build_mapping(permuted_map, oracle_map);
permuted_raw_A = raw_A(permutation, permutation);
declared_A = map_matrix(permuted_raw_A, declared_mapping, pressure, true);
declared_control = compare_matrices(declared_A, mapped_A);

undeclared_A = map_matrix(permuted_raw_A, mapping, pressure, true);
undeclared_permutation_detected = compare_matrices(undeclared_A, oracle_A).e_inf > target;
velocity_sign_A = mapped_A;
velocity_sign_A(~pressure, :) = -velocity_sign_A(~pressure, :);
velocity_sign_detected = compare_matrices(velocity_sign_A, oracle_A).e_inf > target;
missing_pressure_sign_A = map_matrix(raw_A, mapping, pressure, false);
pressure_sign_detected = compare_matrices(missing_pressure_sign_A, oracle_A).e_inf > target;

gauged_state = state;
gauged_state(pressure) = gauged_state(pressure)+3.25;
legal_pressure_gauge_accepted = compare_vectors_gauged(state, gauged_state, pressure).e_inf <= target;
velocity_shifted = state;
velocity_shifted(find(~pressure, 1)) = velocity_shifted(find(~pressure, 1))+1;
illegal_velocity_shift_detected = compare_vectors_gauged(state, velocity_shifted, pressure).e_inf > target;

[~, largest_index] = max(abs(nonzeros(mapped_A)));
[rows, columns, values] = find(mapped_A);
omitted_A = mapped_A;
omitted_A(rows(largest_index), columns(largest_index)) = ...
    omitted_A(rows(largest_index), columns(largest_index))-values(largest_index);
matrix_omission_detected = compare_matrices(omitted_A, oracle_A).e_inf > target;

omitted_patches = patches;
omitted_patches{1} = omitted_patches{1}(2:end);
patch_omission_detected = ~compare_patch_data( ...
    zeros(numel(patches), 1), omitted_patches, zeros(numel(oracle_patches), 1), oracle_patches).membership_and_order_exact;
reordered_patches = patches(end:-1:1);
patch_reorder_detected = ~all(cellfun(@isequal, reordered_patches, oracle_patches));

result = struct( ...
    'declared_permutation_control', declared_report.pass && declared_control.e_inf <= target, ...
    'undeclared_permutation_detected', undeclared_permutation_detected, ...
    'velocity_sign_detected', velocity_sign_detected, ...
    'pressure_equation_sign_detected', pressure_sign_detected, ...
    'legal_pressure_gauge_accepted', legal_pressure_gauge_accepted, ...
    'illegal_velocity_shift_detected', illegal_velocity_shift_detected, ...
    'matrix_entry_omission_detected', matrix_omission_detected, ...
    'patch_dof_omission_detected', patch_omission_detected, ...
    'patch_ordinal_reorder_detected', patch_reorder_detected);
fields = fieldnames(result);
result.pass = all(cellfun(@(name) logical(result.(name)), fields));
end

function output = permute_dof_map(input, permutation)
output = input;
fields = {'kind', 'axis', 'i', 'j', 'position'};
for k = 1:numel(fields)
    name = fields{k};
    output.(name) = input.(name)(permutation, :);
end
end

function result = compare_matrices(left, right)
if ~isequal(size(left), size(right))
    error('Matrix dimensions differ.');
end
difference = left-right;
maximum = full(max(abs(difference(:))));
scale = max([1, full(max(abs(left(:)))), full(max(abs(right(:))))]);
left_pattern = spones(left);
right_pattern = spones(right);
result = struct('max_absolute_error', maximum, 'denominator', scale, ...
                'e_inf', maximum/scale, ...
                'structure_exact', isequal(left_pattern, right_pattern), ...
                'left_nnz', int32(nnz(left)), 'right_nnz', int32(nnz(right)));
end

function result = compare_vectors(left, right)
if numel(left) ~= numel(right)
    error('Vector dimensions differ.');
end
maximum = max(abs(left(:)-right(:)));
scale = max([1; abs(left(:)); abs(right(:))]);
result = struct('max_absolute_error', maximum, 'denominator', scale, ...
                'e_inf', maximum/scale);
end

function result = compare_vectors_gauged(left, right, pressure)
left = left(:);
right = right(:);
left_mean = mean(left(pressure));
right_mean = mean(right(pressure));
left(pressure) = left(pressure)-left_mean;
right(pressure) = right(pressure)-right_mean;
result = compare_vectors(left, right);
result.left_pressure_mean = left_mean;
result.right_pressure_mean = right_mean;
end

function map = read_dof_map(filename)
fid = fopen(filename, 'r');
if fid < 0
    error('Unable to open %s.', filename);
end
cleanup = onCleanup(@() fclose(fid));
header = textscan(fgetl(fid), '%s %d %d');
if isempty(header{1}) || ~strcmp(header{1}{1}, 'ibamr-cav-global-dof-map-v1')
    error('Unsupported DOF-map schema in %s.', filename);
end
map.dimension = double(header{2});
map.count = double(header{3});
records = textscan(fid, '%d %s %d %d %d %f %f');
if numel(records{1}) ~= map.count || ~isequal(double(records{1}), (0:map.count-1)')
    error('Invalid DOF records in %s.', filename);
end
map.kind = records{2};
map.axis = double(records{3});
map.i = double(records{4});
map.j = double(records{5});
map.position = [records{6}, records{7}];
end

function values = read_index_list(filename)
fid = fopen(filename, 'r');
if fid < 0
    error('Unable to open %s.', filename);
end
cleanup = onCleanup(@() fclose(fid));
header = textscan(fgetl(fid), '%s %d');
if isempty(header{1}) || ~strcmp(header{1}{1}, 'ibamr-cav-index-list-v1')
    error('Unsupported index-list schema in %s.', filename);
end
count = double(header{2});
values = fscanf(fid, '%d');
if numel(values) ~= count
    error('Invalid index-list length in %s.', filename);
end
end

function sets = read_index_sets(filename)
fid = fopen(filename, 'r');
if fid < 0
    error('Unable to open %s.', filename);
end
cleanup = onCleanup(@() fclose(fid));
header = textscan(fgetl(fid), '%s %d');
if isempty(header{1}) || ~strcmp(header{1}{1}, 'ibamr-cav-index-sets-v1')
    error('Unsupported index-set schema in %s.', filename);
end
count = double(header{2});
sets = cell(count, 1);
for k = 1:count
    values = sscanf(fgetl(fid), '%d');
    if numel(values) < 2 || values(1) ~= k-1 || numel(values) ~= values(2)+2
        error('Invalid index-set record %d in %s.', k-1, filename);
    end
    sets{k} = double(values(3:end));
end
end

function value = read_matrix_market(filename)
fid = fopen(filename, 'r');
if fid < 0
    error('Unable to open %s.', filename);
end
cleanup = onCleanup(@() fclose(fid));
header = strtrim(fgetl(fid));
dimensions = next_data_line(fid);
dims = sscanf(dimensions, '%d');
if strcmp(header, '%%MatrixMarket matrix coordinate real general')
    if numel(dims) ~= 3
        error('Invalid coordinate dimensions in %s.', filename);
    end
    entries = textscan(fid, '%d %d %f');
    if numel(entries{1}) ~= dims(3)
        error('Invalid coordinate entry count in %s.', filename);
    end
    value = sparse(double(entries{1}), double(entries{2}), entries{3}, dims(1), dims(2));
elseif strcmp(header, '%%MatrixMarket matrix array real general')
    if numel(dims) ~= 2
        error('Invalid array dimensions in %s.', filename);
    end
    entries = textscan(fid, '%f');
    if numel(entries{1}) ~= prod(dims)
        error('Invalid array entry count in %s.', filename);
    end
    value = reshape(entries{1}, dims(1), dims(2));
else
    error('Unsupported MatrixMarket header in %s: %s', filename, header);
end
end

function value = read_vector_market(filename)
value = read_matrix_market(filename);
if size(value, 2) ~= 1
    error('Expected a vector in %s.', filename);
end
value = full(value);
end

function line = next_data_line(fid)
line = fgetl(fid);
while ischar(line) && (isempty(strtrim(line)) || startsWith(strtrim(line), '%'))
    line = fgetl(fid);
end
if ~ischar(line)
    error('Unexpected end of MatrixMarket file.');
end
end

function write_json(filename, value)
fid = fopen(filename, 'w');
if fid < 0
    error('Unable to open %s for writing.', filename);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(value, PrettyPrint=true));
end
