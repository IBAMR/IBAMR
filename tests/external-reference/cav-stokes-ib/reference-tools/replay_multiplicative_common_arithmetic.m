function replay_multiplicative_common_arithmetic()
% Replay the N0--N2 multiplicative-CAV matrix in one arithmetic runtime.

candidate_root = getenv('CAV_CANDIDATE_ROOT');
oracle_root = getenv('CAV_ORACLE_ROOT');
scope = getenv('CAV_REPLAY_SCOPE');
report_path = getenv('CAV_REPLAY_REPORT');
if isempty(candidate_root) || isempty(oracle_root) || isempty(scope) || isempty(report_path)
    error(['CAV_CANDIDATE_ROOT, CAV_ORACLE_ROOT, CAV_REPLAY_SCOPE, and ' ...
           'CAV_REPLAY_REPORT are required.']);
end
if strcmp(scope, 'N2')
    oracle_worktree = getenv('CAV_ORACLE_WORKTREE');
    if isempty(oracle_worktree)
        error('CAV_ORACLE_WORKTREE is required for N2 solver utilities.');
    end
    addpath(oracle_worktree);
    setup_project_paths(oracle_worktree);
end
if strcmp(scope, 'N0')
    stiffnesses = 1;
    directions = {'forward'};
elseif strcmp(scope, 'N1')
    stiffnesses = [1, 1e2, 1e4];
    directions = {'forward', 'reverse'};
elseif strcmp(scope, 'N2')
    stiffnesses = [1, 1e2, 1e4];
    directions = {'forward'};
else
    error('CAV_REPLAY_SCOPE must be N0, N1, or N2.');
end

cases = struct([]);
for k = 1:numel(stiffnesses)
    if strcmp(scope, 'N2')
        current_case = replay_stiffness_case_n2( ...
            candidate_root, oracle_root, stiffnesses(k));
    else
        current_case = replay_stiffness_case( ...
            candidate_root, oracle_root, stiffnesses(k), directions);
    end
    if isempty(cases)
        cases = current_case;
    else
        cases(k) = current_case;
    end
end

report = struct();
report.schema = 'cav-multiplicative-paired-common-arithmetic-v4';
report.scope = scope;
if strcmp(scope, 'N2')
    report.supported_geometry = ['live periodic 8x8-to-16x16 hierarchy; both legal ' ...
        'SAMRAI levels participate in the two-level replay'];
    report.traversal_contract = ['the normative sandbox and live production V-cycle uses ' ...
        'forward pre- and forward post-sweeps; a forward-pre/reverse-post palindrome is ' ...
        'a separate sensitivity replay and is not a production option'];
else
    report.supported_geometry = ['live periodic 4x4-to-8x8 hierarchy; the finest live ' ...
        '8x8 pressure-CAV level is replayed as the transfer-free N1 level'];
end
report.strict_authority = ['independent IBAMR STRICT extension: retain candidate cells only when ' ...
    'their complete incident-velocity stencil is contained in the E_h-expanded velocity set'];
report.reverse_policy = ['reverse is a sensitivity replay of the exported patch sequence, not a ' ...
    'production option or a second stored patch representation'];
report.error_formula = ['E_inf=max(abs(x-y))/max(1,max(abs(x)),max(abs(y))); ' ...
    'matrices use the same denominator floor and all entries'];
report.runtime = struct('name', 'MATLAB', 'version', version, 'computer', computer);
report.cases = num2cell(cases);
report.pass = all([cases.pass]);

write_json(report_path, report);
for k = 1:numel(cases)
    fprintf('%s_K=%.17g_operator_e_inf=%.17e\n', lower(scope), ...
        cases(k).physical_stiffness_K, cases(k).operator.e_inf);
    fprintf('%s_K=%.17g_relaxed_patch_exact=%d\n', lower(scope), ...
        cases(k).physical_stiffness_K, cases(k).patches.relaxed_candidate_vs_oracle.pass);
    if ~strcmp(scope, 'N2')
        fprintf('%s_K=%.17g_strict_patch_exact=%d\n', lower(scope), ...
            cases(k).physical_stiffness_K, cases(k).patches.strict_candidate_vs_independent.pass);
    else
        fprintf('%s_K=%.17g_transfer_actions_pass=%d\n', lower(scope), ...
            cases(k).physical_stiffness_K, cases(k).live_fac_stage_validation.pass);
        fprintf('%s_K=%.17g_vcycle_pass=%d\n', lower(scope), ...
            cases(k).physical_stiffness_K, cases(k).normative_forward_forward.pass);
        fprintf('%s_K=%.17g_fgmres_pass=%d\n', lower(scope), ...
            cases(k).physical_stiffness_K, cases(k).common_arithmetic_fgmres.pass);
    end
end
fprintf('%s_paired_replay_pass=%d\n', lower(scope), report.pass);
if ~report.pass
    error('%s paired common-arithmetic replay failed; see %s.', scope, report_path);
end
end

function result = replay_stiffness_case(candidate_root, oracle_root, stiffness, directions)
label = stiffness_label(stiffness);
target = fixed_target(stiffness);
relaxed_dir = fullfile(candidate_root, label, 'RELAXED');
strict_dir = fullfile(candidate_root, label, 'STRICT');
oracle_dir = fullfile(oracle_root, label);

relaxed_map = read_dof_map(fullfile(relaxed_dir, 'cav_live_export_contract_dof_map.txt'));
strict_map = read_dof_map(fullfile(strict_dir, 'cav_live_export_contract_dof_map.txt'));
oracle_map = read_dof_map(fullfile(oracle_dir, 'oracle_dof_map.txt'));
[mapping, mapping_report] = build_mapping(relaxed_map, oracle_map);
[strict_mapping, strict_mapping_report] = build_mapping(strict_map, oracle_map);
mapping_report.strict_map_same = strict_mapping_report.pass && isequal(mapping, strict_mapping);
mapping_report.pass = mapping_report.pass && mapping_report.strict_map_same;
pressure = strcmp(oracle_map.kind, 'P');
configuration = validate_configuration( ...
    relaxed_dir, strict_dir, oracle_dir, stiffness, relaxed_map, oracle_map);

relaxed_A_raw = read_matrix_market(fullfile(relaxed_dir, 'cav_live_export_contract_A.mtx'));
strict_A_raw = read_matrix_market(fullfile(strict_dir, 'cav_live_export_contract_A.mtx'));
relaxed_E_raw = read_matrix_market(fullfile(relaxed_dir, 'cav_live_export_contract_E_h.mtx'));
strict_E_raw = read_matrix_market(fullfile(strict_dir, 'cav_live_export_contract_E_h.mtx'));
candidate_A = map_matrix(relaxed_A_raw, mapping, pressure, true);
candidate_E = map_matrix(relaxed_E_raw, mapping, pressure, false);
strict_A = map_matrix(strict_A_raw, strict_mapping, pressure, true);
strict_E = map_matrix(strict_E_raw, strict_mapping, pressure, false);
oracle_A = read_matrix_market(fullfile(oracle_dir, 'oracle_A.mtx'));
oracle_E = read_matrix_market(fullfile(oracle_dir, 'oracle_E_h.mtx'));

operator_comparison = compare_matrices(candidate_A, oracle_A);
elasticity_comparison = compare_matrices(candidate_E, oracle_E);
policy_operator_control = compare_matrices(candidate_A, strict_A);
policy_elasticity_control = compare_matrices(candidate_E, strict_E);

relaxed_seeds_raw = read_index_list(fullfile(relaxed_dir, 'cav_live_export_contract_pressure_seeds.txt'));
strict_seeds_raw = read_index_list(fullfile(strict_dir, 'cav_live_export_contract_pressure_seeds.txt'));
relaxed_patches_raw = read_index_sets(fullfile(relaxed_dir, 'cav_live_export_contract_patches.txt'));
strict_patches_raw = read_index_sets(fullfile(strict_dir, 'cav_live_export_contract_patches.txt'));
relaxed_seeds = mapping(relaxed_seeds_raw+1);
strict_seeds = strict_mapping(strict_seeds_raw+1);
relaxed_patches = map_sets(relaxed_patches_raw, mapping);
strict_patches = map_sets(strict_patches_raw, strict_mapping);
oracle_seeds = read_index_list(fullfile(oracle_dir, 'oracle_pressure_seeds.txt'))+1;
oracle_relaxed_patches = one_based_sets(read_index_sets(fullfile(oracle_dir, 'oracle_patches.txt')));
independent_relaxed_patches = construct_reference_patches( ...
    oracle_E, oracle_map, oracle_seeds, 'RELAXED');
independent_strict_patches = construct_reference_patches( ...
    oracle_E, oracle_map, oracle_seeds, 'STRICT');

patches = struct();
patches.relaxed_candidate_vs_oracle = compare_patch_data( ...
    relaxed_seeds, relaxed_patches, oracle_seeds, oracle_relaxed_patches);
patches.relaxed_oracle_vs_independent = compare_patch_data( ...
    oracle_seeds, oracle_relaxed_patches, oracle_seeds, independent_relaxed_patches);
patches.strict_candidate_vs_independent = compare_patch_data( ...
    strict_seeds, strict_patches, oracle_seeds, independent_strict_patches);

candidate_rhs_relaxed = read_candidate_state(relaxed_dir, mapping, pressure, ...
    'cav_fac_stage0_pre-smooth-input_level1_rhs.mtx', true);
candidate_initial_relaxed = read_candidate_state(relaxed_dir, mapping, pressure, ...
    'cav_fac_stage0_pre-smooth-input_level1_solution.mtx', false);
candidate_actual_relaxed = read_candidate_state(relaxed_dir, mapping, pressure, ...
    'cav_fac_stage1_pre-smooth-output_level1_solution.mtx', false);
candidate_rhs_strict = read_candidate_state(strict_dir, strict_mapping, pressure, ...
    'cav_fac_stage0_pre-smooth-input_level1_rhs.mtx', true);
candidate_initial_strict = read_candidate_state(strict_dir, strict_mapping, pressure, ...
    'cav_fac_stage0_pre-smooth-input_level1_solution.mtx', false);
candidate_actual_strict = read_candidate_state(strict_dir, strict_mapping, pressure, ...
    'cav_fac_stage1_pre-smooth-output_level1_solution.mtx', false);
oracle_rhs = read_vector_market(fullfile(oracle_dir, 'oracle_rhs.mtx'));

policy_input_control = struct( ...
    'rhs', compare_vectors(candidate_rhs_relaxed, candidate_rhs_strict), ...
    'initial', compare_vectors(candidate_initial_relaxed, candidate_initial_strict));
policy_input_control.pass = policy_input_control.rhs.e_inf == 0 && ...
                            policy_input_control.initial.e_inf == 0;

relaxed = replay_policy(candidate_A, candidate_rhs_relaxed, candidate_initial_relaxed, ...
    candidate_actual_relaxed, relaxed_patches, oracle_relaxed_patches, oracle_A, oracle_rhs, ...
    pressure, target, directions);
strict = replay_policy(candidate_A, candidate_rhs_strict, candidate_initial_strict, ...
    candidate_actual_strict, strict_patches, independent_strict_patches, oracle_A, oracle_rhs, ...
    pressure, target, directions);
relaxed.native_all_patch_diagnostics = read_native_diagnostics(relaxed_dir, numel(relaxed_patches));
strict.native_all_patch_diagnostics = read_native_diagnostics(strict_dir, numel(strict_patches));
relaxed.candidate_selected_local_solve = check_candidate_local_trace( ...
    relaxed_dir, relaxed_A_raw, relaxed_patches_raw{1}, target);
strict.candidate_selected_local_solve = check_candidate_local_trace( ...
    strict_dir, strict_A_raw, strict_patches_raw{1}, target);

[~, ~, oracle_patch0] = apply_sweep( ...
    oracle_A, oracle_rhs, zeros(size(oracle_rhs)), oracle_relaxed_patches, pressure);
[oracle_replay_w, oracle_replay_r] = apply_sweep( ...
    oracle_A, oracle_rhs, zeros(size(oracle_rhs)), oracle_relaxed_patches, pressure);
oracle_actual = read_vector_market(fullfile(oracle_dir, 'oracle_smoother_correction.mtx'));
oracle_actual_r = read_vector_market(fullfile(oracle_dir, 'oracle_smoother_fresh_residual.mtx'));
oracle_pinned_control = struct( ...
    'correction', compare_vectors_gauged(oracle_replay_w, oracle_actual, pressure), ...
    'fresh_residual', compare_vectors(oracle_replay_r, oracle_actual_r));
oracle_pinned_control.within_fixed_target = ...
    oracle_pinned_control.correction.e_inf <= target && ...
    oracle_pinned_control.fresh_residual.e_inf <= target;
oracle_pinned_control.classification = ...
    'reported cross-runtime backend/arithmetic diagnostic; not a common-arithmetic gate';
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

mutations = run_mutations(relaxed_A_raw, candidate_A, oracle_A, relaxed_map, ...
    oracle_map, mapping, pressure, relaxed_patches, oracle_relaxed_patches, ...
    relaxed.forward_state, target);

result = struct();
result.physical_stiffness_K = stiffness;
result.fixed_target = target;
result.configuration = configuration;
result.mapping = mapping_report;
result.operator = operator_comparison;
result.elasticity = elasticity_comparison;
result.policy_operator_identity = policy_operator_control;
result.policy_elasticity_identity = policy_elasticity_control;
result.policy_input_identity = policy_input_control;
result.patches = patches;
result.relaxed = rmfield(relaxed, 'forward_state');
result.strict = rmfield(strict, 'forward_state');
result.pinned_oracle_forward_diagnostic_not_gate = oracle_pinned_control;
result.oracle_selected_local_solve = oracle_local_control;
result.mutations = mutations;
result.pass = configuration.pass && mapping_report.pass && matrix_within(operator_comparison, target) && ...
    matrix_within(elasticity_comparison, target) && matrix_exact(policy_operator_control) && ...
    matrix_exact(policy_elasticity_control) && policy_input_control.pass && ...
    patches.relaxed_candidate_vs_oracle.pass && patches.relaxed_oracle_vs_independent.pass && ...
    patches.strict_candidate_vs_independent.pass && relaxed.pass && strict.pass && ...
    relaxed.native_all_patch_diagnostics.pass && strict.native_all_patch_diagnostics.pass && ...
    relaxed.candidate_selected_local_solve.pass && strict.candidate_selected_local_solve.pass && ...
    oracle_local_control.pass && mutations.pass;
end

function result = replay_stiffness_case_n2(candidate_root, oracle_root, stiffness)
label = stiffness_label(stiffness);
target = fixed_target(stiffness);
candidate_dir = fullfile(candidate_root, label, 'RELAXED');
oracle_dir = fullfile(oracle_root, label);

candidate_map = read_dof_map(fullfile(candidate_dir, 'cav_live_export_contract_dof_map.txt'));
oracle_map = read_dof_map(fullfile(oracle_dir, 'oracle_dof_map.txt'));
[mapping, mapping_report] = build_mapping(candidate_map, oracle_map);
pressure = strcmp(oracle_map.kind, 'P');
candidate_coarse_map = read_dof_map(fullfile(candidate_dir, 'cav_dynamic_level0_dof_map.txt'));
oracle_coarse_map = read_dof_map(fullfile(oracle_dir, 'oracle_coarse_dof_map.txt'));
[coarse_mapping, coarse_mapping_report] = build_mapping(candidate_coarse_map, oracle_coarse_map);
coarse_pressure = strcmp(oracle_coarse_map.kind, 'P');

candidate_A_raw = read_matrix_market(fullfile(candidate_dir, 'cav_live_export_contract_A.mtx'));
candidate_E_raw = read_matrix_market(fullfile(candidate_dir, 'cav_live_export_contract_E_h.mtx'));
candidate_coarse_A_raw = read_matrix_market(fullfile(candidate_dir, 'cav_dynamic_level0_A.mtx'));
candidate_A = map_matrix(candidate_A_raw, mapping, pressure, true);
candidate_E = map_matrix(candidate_E_raw, mapping, pressure, false);
candidate_coarse_A = map_matrix( ...
    candidate_coarse_A_raw, coarse_mapping, coarse_pressure, true);
oracle_A = read_matrix_market(fullfile(oracle_dir, 'oracle_A.mtx'));
oracle_E = read_matrix_market(fullfile(oracle_dir, 'oracle_E_h.mtx'));
oracle_coarse_A = read_matrix_market(fullfile(oracle_dir, 'oracle_coarse_A.mtx'));
restriction = read_matrix_market(fullfile(oracle_dir, 'oracle_restriction.mtx'));
prolongation = read_matrix_market(fullfile(oracle_dir, 'oracle_prolongation.mtx'));

operator_comparison = compare_matrices(candidate_A, oracle_A);
elasticity_comparison = compare_matrices(candidate_E, oracle_E);
coarse_operator_comparison = compare_matrices(candidate_coarse_A, oracle_coarse_A);

candidate_seeds_raw = read_index_list(fullfile( ...
    candidate_dir, 'cav_live_export_contract_pressure_seeds.txt'));
candidate_patches_raw = read_index_sets(fullfile( ...
    candidate_dir, 'cav_live_export_contract_patches.txt'));
candidate_seeds = mapping(candidate_seeds_raw+1);
candidate_patches = map_sets(candidate_patches_raw, mapping);
oracle_seeds = read_index_list(fullfile(oracle_dir, 'oracle_pressure_seeds.txt'))+1;
oracle_patches = one_based_sets(read_index_sets(fullfile(oracle_dir, 'oracle_patches.txt')));
independent_patches = construct_reference_patches( ...
    oracle_E, oracle_map, oracle_seeds, 'RELAXED');
patches = struct( ...
    'relaxed_candidate_vs_oracle', compare_patch_data( ...
        candidate_seeds, candidate_patches, oracle_seeds, oracle_patches), ...
    'relaxed_oracle_vs_independent', compare_patch_data( ...
        oracle_seeds, oracle_patches, oracle_seeds, independent_patches));

configuration = validate_configuration_n2( ...
    candidate_dir, oracle_dir, stiffness, candidate_map, oracle_map, ...
    candidate_coarse_map, oracle_coarse_map);

fine_rhs = read_candidate_state(candidate_dir, mapping, pressure, ...
    'cav_fac_stage0_pre-smooth-input_level1_rhs.mtx', true);
pre_input = read_candidate_state(candidate_dir, mapping, pressure, ...
    'cav_fac_stage0_pre-smooth-input_level1_solution.mtx', false);
pre_output = read_candidate_state(candidate_dir, mapping, pressure, ...
    'cav_fac_stage1_pre-smooth-output_level1_solution.mtx', false);
coarse_rhs = read_candidate_state(candidate_dir, coarse_mapping, coarse_pressure, ...
    'cav_fac_stage2_coarse-rhs_level0_rhs.mtx', true);
coarse_input = read_candidate_state(candidate_dir, coarse_mapping, coarse_pressure, ...
    'cav_fac_stage2_coarse-rhs_level0_solution.mtx', false);
coarse_output = read_candidate_state(candidate_dir, coarse_mapping, coarse_pressure, ...
    'cav_fac_stage3_coarse-correction_level0_solution.mtx', false);
post_input = read_candidate_state(candidate_dir, mapping, pressure, ...
    'cav_fac_stage4_post-smooth-input_level1_solution.mtx', false);
post_output = read_candidate_state(candidate_dir, mapping, pressure, ...
    'cav_fac_stage5_post-smooth-output_level1_solution.mtx', false);
post_rhs = read_candidate_state(candidate_dir, mapping, pressure, ...
    'cav_fac_stage5_post-smooth-output_level1_rhs.mtx', true);

[replayed_pre, replayed_pre_residual] = apply_ordered_sweep( ...
    candidate_A, fine_rhs, pre_input, candidate_patches, pressure, 'forward');
replayed_coarse_rhs = restriction*replayed_pre_residual;
replayed_coarse = solve_mean_zero_pressure( ...
    candidate_coarse_A, replayed_coarse_rhs, coarse_pressure);
replayed_post_input = replayed_pre+prolongation*replayed_coarse;
[replayed_post, replayed_post_residual] = apply_ordered_sweep( ...
    candidate_A, fine_rhs, replayed_post_input, candidate_patches, pressure, 'forward');

live_fac_stage_validation = struct( ...
    'fine_rhs_unchanged', compare_vectors(fine_rhs, post_rhs), ...
    'pre_input_zero', compare_vectors(pre_input, zeros(size(pre_input))), ...
    'pre_output', compare_vectors_gauged(replayed_pre, pre_output, pressure), ...
    'restricted_coarse_rhs', compare_vectors(replayed_coarse_rhs, coarse_rhs), ...
    'coarse_input_zero', compare_vectors(coarse_input, zeros(size(coarse_input))), ...
    'coarse_correction', compare_vectors_gauged(replayed_coarse, coarse_output, coarse_pressure), ...
    'prolonged_post_input', compare_vectors_gauged(replayed_post_input, post_input, pressure), ...
    'post_output', compare_vectors_gauged(replayed_post, post_output, pressure), ...
    'post_fresh_residual', compare_vectors(replayed_post_residual, fine_rhs-candidate_A*post_output));
stage_fields = fieldnames(live_fac_stage_validation);
stage_errors = cellfun(@(name) live_fac_stage_validation.(name).e_inf, stage_fields);
live_fac_stage_validation.within_fixed_target = all(stage_errors <= target);
live_fac_stage_validation.within_hard_ceiling = all(stage_errors <= 1e-10);
live_fac_stage_validation.classification = [ ...
    'live-versus-common-runtime stage differences are retained diagnostics; ' ...
    'fixed per-K targets apply to the paired common-arithmetic result'];
live_fac_stage_validation.pass = ...
    live_fac_stage_validation.fine_rhs_unchanged.e_inf == 0 && ...
    live_fac_stage_validation.pre_input_zero.e_inf == 0 && ...
    live_fac_stage_validation.coarse_input_zero.e_inf == 0 && ...
    live_fac_stage_validation.within_hard_ceiling;

restriction_mutated = restriction;
[r_row, r_column] = find(restriction_mutated, 1);
restriction_mutated(r_row, r_column) = -restriction_mutated(r_row, r_column);
prolongation_mutated = prolongation;
[p_row, p_column] = find(prolongation_mutated, 1);
prolongation_mutated(p_row, p_column) = -prolongation_mutated(p_row, p_column);
transfer_mutations = struct( ...
    'restriction_sign_detected', compare_vectors( ...
        restriction_mutated*replayed_pre_residual, coarse_rhs).e_inf > target, ...
    'prolongation_sign_detected', compare_vectors_gauged( ...
        replayed_pre+prolongation_mutated*replayed_coarse, post_input, pressure).e_inf > target);
transfer_mutations.pass = transfer_mutations.restriction_sign_detected && ...
                          transfer_mutations.prolongation_sign_detected;

oracle_rhs = read_vector_market(fullfile(oracle_dir, 'oracle_rhs.mtx'));
normative_forward_forward = paired_vcycle_on_both_inputs( ...
    candidate_A, candidate_coarse_A, fine_rhs, candidate_patches, oracle_patches, ...
    oracle_A, oracle_coarse_A, oracle_rhs, restriction, prolongation, ...
    pressure, coarse_pressure, target, 'forward');
palindromic_sensitivity = paired_vcycle_on_both_inputs( ...
    candidate_A, candidate_coarse_A, fine_rhs, candidate_patches, oracle_patches, ...
    oracle_A, oracle_coarse_A, oracle_rhs, restriction, prolongation, ...
    pressure, coarse_pressure, target, 'reverse');
palindromic_difference = paired_result( ...
    normative_forward_forward.live_candidate_state, ...
    normative_forward_forward.live_candidate_residual, ...
    palindromic_sensitivity.live_candidate_state, ...
    palindromic_sensitivity.live_candidate_residual, pressure, target);
palindromic_sensitivity.differs_from_normative = ...
    palindromic_difference.correction.e_inf > target || ...
    palindromic_difference.fresh_residual.e_inf > target;
palindromic_sensitivity.normative_comparison = palindromic_difference;
palindromic_sensitivity.pass = palindromic_sensitivity.pass && ...
                               palindromic_sensitivity.differs_from_normative;
normative_forward_forward = rmfield(normative_forward_forward, ...
    {'live_candidate_state', 'live_candidate_residual'});
palindromic_sensitivity = rmfield(palindromic_sensitivity, ...
    {'live_candidate_state', 'live_candidate_residual'});

common_arithmetic_fgmres = paired_fgmres_on_both_inputs( ...
    candidate_A, candidate_coarse_A, fine_rhs, candidate_patches, oracle_patches, ...
    oracle_A, oracle_coarse_A, oracle_rhs, restriction, prolongation, ...
    pressure, coarse_pressure, target);

live_production = validate_live_production_n2( ...
    candidate_dir, oracle_dir, mapping, pressure, candidate_A, fine_rhs, ...
    common_arithmetic_fgmres.live_candidate_solution);
common_arithmetic_fgmres = rmfield(common_arithmetic_fgmres, ...
    {'live_candidate_solution'});

native_diagnostics = read_native_diagnostics(candidate_dir, numel(candidate_patches));
selected_local_solve = check_candidate_local_trace( ...
    candidate_dir, candidate_A_raw, candidate_patches_raw{1}, target);
mutations = run_mutations(candidate_A_raw, candidate_A, oracle_A, candidate_map, ...
    oracle_map, mapping, pressure, candidate_patches, oracle_patches, replayed_pre, target);

result = struct();
result.physical_stiffness_K = stiffness;
result.fixed_target = target;
result.configuration = configuration;
result.mapping = mapping_report;
result.coarse_mapping = coarse_mapping_report;
result.operator = operator_comparison;
result.elasticity = elasticity_comparison;
result.coarse_operator = coarse_operator_comparison;
result.patches = patches;
result.live_fac_stage_validation = live_fac_stage_validation;
result.transfer_mutations = transfer_mutations;
result.normative_forward_forward = normative_forward_forward;
result.palindromic_forward_reverse_sensitivity = palindromic_sensitivity;
result.common_arithmetic_fgmres = common_arithmetic_fgmres;
result.live_production = live_production;
result.native_all_patch_diagnostics = native_diagnostics;
result.candidate_selected_local_solve = selected_local_solve;
result.mutations = mutations;
result.pass = configuration.pass && mapping_report.pass && coarse_mapping_report.pass && ...
    matrix_within(operator_comparison, target) && matrix_within(elasticity_comparison, target) && ...
    matrix_within(coarse_operator_comparison, target) && ...
    patches.relaxed_candidate_vs_oracle.pass && patches.relaxed_oracle_vs_independent.pass && ...
    live_fac_stage_validation.pass && transfer_mutations.pass && ...
    normative_forward_forward.pass && palindromic_sensitivity.pass && ...
    common_arithmetic_fgmres.pass && live_production.pass && ...
    native_diagnostics.pass && selected_local_solve.pass && mutations.pass;
end

function result = validate_configuration_n2(candidate_dir, oracle_dir, stiffness, ...
                                            candidate_map, oracle_map, ...
                                            candidate_coarse_map, oracle_coarse_map)
manifest = read_live_manifest(fullfile( ...
    candidate_dir, 'cav_live_export_contract_manifest.txt'));
metadata = jsondecode(fileread(fullfile(oracle_dir, 'oracle_metadata.json')));
stdout = fileread(fullfile(candidate_dir, 'candidate.stdout'));
result = struct( ...
    'dimension_2d', candidate_map.dimension == 2 && oracle_map.dimension == 2 && ...
        candidate_coarse_map.dimension == 2 && oracle_coarse_map.dimension == 2 && metadata.dimension == 2, ...
    'fine_global_dof_count_exact', candidate_map.count == 768 && oracle_map.count == 768, ...
    'coarse_global_dof_count_exact', candidate_coarse_map.count == 192 && oracle_coarse_map.count == 192, ...
    'live_hierarchy_is_legal_8_to_16', metadata.base_grid_cells_per_axis == 8 && ...
        metadata.refinement_ratio == 2 && metadata.finest_grid_cells_per_axis == 16, ...
    'single_rank', str2double(manifest.mpi_ranks) == 1, ...
    'physical_stiffness_exact', read_output_number(stdout, 'physical_stiffness') == stiffness && ...
        metadata.physical_stiffness_K == stiffness, ...
    'relaxed_policy_exact', strcmp(manifest.closure_policy, 'RELAXED'), ...
    'seed_and_order_exact', strcmp(manifest.patch_seed_type, 'PRESSURE_CELL') && ...
        strcmp(manifest.seed_stride, '1') && strcmp(manifest.traversal_order, 'I_J'), ...
    'composition_and_backend_exact', strcmp(manifest.composition, 'multiplicative') && ...
        strcmp(manifest.local_solver_backend, 'blas-lapack-lu'), ...
    'sign_and_gauge_exact', strcmp(manifest.pressure_equation, 'minus-div') && ...
        str2double(manifest.pressure_equation_row_multiplier_to_oracle) == -1 && ...
        strcmp(manifest.pressure_gauge, 'zero-mean-correction') && ...
        strcmp(metadata.pressure_equation, 'div') && strcmp(metadata.pressure_gauge, 'zero-mean correction'), ...
    'oracle_contract_exact', strcmp(metadata.oracle_sha, ...
        '5b77344db6746269f8c77695c99e9043907ba74b') && ...
        strcmp(metadata.patch_strategy, 'targeted_ib') && ...
        strcmp(metadata.composition, 'multiplicative') && ...
        strcmp(metadata.delta_function, 'IB_4') && metadata.patch_count == 256);
fields = fieldnames(result);
result.pass = all(cellfun(@(name) logical(result.(name)), fields));
end

function result = paired_vcycle_on_both_inputs(live_A, live_coarse_A, live_rhs, ...
        candidate_patches, reference_patches, oracle_A, oracle_coarse_A, oracle_rhs, ...
        restriction, prolongation, pressure, coarse_pressure, target, post_direction)
[live_candidate, live_candidate_r] = apply_common_vcycle( ...
    live_A, live_coarse_A, live_rhs, candidate_patches, restriction, prolongation, ...
    pressure, coarse_pressure, 'forward', post_direction);
[live_reference, live_reference_r] = apply_common_vcycle( ...
    live_A, live_coarse_A, live_rhs, reference_patches, restriction, prolongation, ...
    pressure, coarse_pressure, 'forward', post_direction);
[oracle_candidate, oracle_candidate_r] = apply_common_vcycle( ...
    oracle_A, oracle_coarse_A, oracle_rhs, candidate_patches, restriction, prolongation, ...
    pressure, coarse_pressure, 'forward', post_direction);
[oracle_reference, oracle_reference_r] = apply_common_vcycle( ...
    oracle_A, oracle_coarse_A, oracle_rhs, reference_patches, restriction, prolongation, ...
    pressure, coarse_pressure, 'forward', post_direction);
result = struct( ...
    'pre_direction', 'forward', 'post_direction', post_direction, ...
    'live_operator', paired_result(live_candidate, live_candidate_r, ...
        live_reference, live_reference_r, pressure, target), ...
    'oracle_operator', paired_result(oracle_candidate, oracle_candidate_r, ...
        oracle_reference, oracle_reference_r, pressure, target), ...
    'live_candidate_state', live_candidate, ...
    'live_candidate_residual', live_candidate_r);
result.pass = result.live_operator.pass && result.oracle_operator.pass;
end

function [w, fresh_residual] = apply_common_vcycle(A, coarse_A, b, patches, ...
        restriction, prolongation, pressure, coarse_pressure, pre_direction, post_direction)
[w, residual] = apply_ordered_sweep( ...
    A, b, zeros(size(b)), patches, pressure, pre_direction);
coarse_rhs = restriction*residual;
coarse_correction = solve_mean_zero_pressure(coarse_A, coarse_rhs, coarse_pressure);
w = w+prolongation*coarse_correction;
[w, fresh_residual] = apply_ordered_sweep( ...
    A, b, w, patches, pressure, post_direction);
end

function [w, fresh_residual] = apply_ordered_sweep(A, b, w, patches, pressure, direction)
if strcmp(direction, 'reverse')
    patches = patches(end:-1:1);
elseif ~strcmp(direction, 'forward')
    error('Unsupported sweep direction %s.', direction);
end
[w, fresh_residual] = apply_sweep(A, b, w, patches, pressure);
end

function correction = solve_mean_zero_pressure(A, rhs, pressure)
constraint = zeros(size(rhs));
constraint(pressure) = 1/nnz(pressure);
augmented = [full(A), constraint; constraint', 0];
solution = augmented\[rhs; 0];
correction = solution(1:numel(rhs));
correction(pressure) = correction(pressure)-mean(correction(pressure));
end

function result = paired_fgmres_on_both_inputs(live_A, live_coarse_A, live_rhs, ...
        candidate_patches, reference_patches, oracle_A, oracle_coarse_A, oracle_rhs, ...
        restriction, prolongation, pressure, coarse_pressure, target)
[live_result, live_candidate_solution] = paired_fgmres( ...
    live_A, live_coarse_A, live_rhs, candidate_patches, reference_patches, ...
    restriction, prolongation, pressure, coarse_pressure, target);
[oracle_result, ~] = paired_fgmres( ...
    oracle_A, oracle_coarse_A, oracle_rhs, candidate_patches, reference_patches, ...
    restriction, prolongation, pressure, coarse_pressure, target);
result = struct('live_operator', live_result, 'oracle_operator', oracle_result, ...
                'live_candidate_solution', live_candidate_solution);
result.pass = live_result.pass && oracle_result.pass;
end

function [result, candidate_solution] = paired_fgmres(A, coarse_A, rhs, ...
        candidate_patches, reference_patches, restriction, prolongation, ...
        pressure, coarse_pressure, target)
candidate_pc = @(input) apply_common_vcycle(A, coarse_A, input, candidate_patches, ...
    restriction, prolongation, pressure, coarse_pressure, 'forward', 'forward');
reference_pc = @(input) apply_common_vcycle(A, coarse_A, input, reference_patches, ...
    restriction, prolongation, pressure, coarse_pressure, 'forward', 'forward');
[candidate_solution, candidate_flag, candidate_relres, candidate_iterations, candidate_history] = ...
    fgmres_opt(A, rhs, candidate_pc, zeros(size(rhs)), 1e-10, 1e-50, 200, false, 'either');
[reference_solution, reference_flag, reference_relres, reference_iterations, reference_history] = ...
    fgmres_opt(A, rhs, reference_pc, zeros(size(rhs)), 1e-10, 1e-50, 200, false, 'either');
candidate_residual = rhs-A*candidate_solution;
reference_residual = rhs-A*reference_solution;
result = struct( ...
    'solution', compare_vectors_gauged(candidate_solution, reference_solution, pressure), ...
    'fresh_residual', compare_vectors(candidate_residual, reference_residual), ...
    'relative_history', compare_histories(candidate_history, reference_history), ...
    'candidate_flag', int32(candidate_flag), 'reference_flag', int32(reference_flag), ...
    'candidate_iterations', int32(candidate_iterations), ...
    'reference_iterations', int32(reference_iterations), ...
    'candidate_reported_relative_residual', candidate_relres, ...
    'reference_reported_relative_residual', reference_relres);
result.pass = candidate_flag == reference_flag && candidate_iterations == reference_iterations && ...
    result.solution.e_inf <= target && result.fresh_residual.e_inf <= target && ...
    result.relative_history.e_inf <= target;
end

function result = validate_live_production_n2(candidate_dir, oracle_dir, mapping, pressure, ...
                                              A, rhs, common_solution)
trace = read_krylov_trace(fullfile(candidate_dir, 'cav_fgmres_trace.txt'));
last_solution = read_candidate_state(candidate_dir, mapping, pressure, ...
    [trace.last_stem '_level1_solution.mtx'], false);
candidate_physical = read_named_scalars(fullfile(candidate_dir, 'cav_physical_residual.txt'));
oracle_summary = jsondecode(fileread(fullfile(oracle_dir, 'oracle_fgmres_summary.json')));
oracle_solution = read_vector_market(fullfile(oracle_dir, 'oracle_fgmres_solution.mtx'));
oracle_rhs = read_vector_market(fullfile(oracle_dir, 'oracle_rhs.mtx'));
oracle_history = double(oracle_summary.residual_history(:));
candidate_scalar = [candidate_physical.relative_residual; ...
                    candidate_physical.residual_momentum_l2; ...
                    candidate_physical.residual_divergence_l2];
oracle_scalar = [double(oracle_summary.relative_residual); ...
                 double(oracle_summary.residual_momentum_l2); ...
                 double(oracle_summary.residual_divergence_l2)];
candidate_reporting_values = [candidate_scalar; ...
    candidate_physical.ib_matrix_free_action_l2; ...
    candidate_physical.ib_matrix_action_l2; ...
    candidate_physical.ib_residual_l2; ...
    candidate_physical.ib_relative_residual];
candidate_reporting_complete = all(isfinite(candidate_reporting_values));
oracle_reporting_complete = all(isfinite([oracle_scalar; oracle_history]));
history_reporting_complete = ~isempty(trace.residuals) && ...
    all(isfinite(trace.residuals)) && ~isempty(oracle_history);
candidate_original_residual_pass = trace.count >= 2 && candidate_reporting_complete && ...
    candidate_physical.relative_residual <= 1e-10;
result = struct( ...
    'candidate_converged', candidate_original_residual_pass, ...
    'candidate_original_residual_target', 1e-10, ...
    'candidate_original_residual_pass', candidate_original_residual_pass, ...
    'candidate_physical_reporting_complete', candidate_reporting_complete, ...
    'oracle_physical_reporting_complete', oracle_reporting_complete, ...
    'iteration_history_reporting_complete', history_reporting_complete, ...
    'candidate_ib_consistency_previous_fixed_target_diagnostic', ...
        candidate_physical.ib_relative_residual <= 1e-10, ...
    'candidate_iteration_count', int32(trace.count-1), ...
    'oracle_iteration_count', int32(oracle_summary.iterations), ...
    'iteration_count_exact', trace.count-1 == double(oracle_summary.iterations), ...
    'normalized_history', compare_histories(trace.residuals, oracle_history), ...
    'physical_scalar', compare_vectors(candidate_scalar, oracle_scalar), ...
    'solution', compare_vectors_gauged(last_solution, oracle_solution, pressure), ...
    'rhs', compare_vectors(rhs, oracle_rhs), ...
    'common_replay_vs_live_solution', compare_vectors_gauged(common_solution, last_solution, pressure), ...
    'common_replay_fresh_residual_vs_live', compare_vectors( ...
        rhs-A*common_solution, rhs-A*last_solution));
result.cross_runtime_classification = [ ...
    'candidate/oracle physical scalars, iteration count/history, rho_IB, raw RHS, final state, ' ...
    'and common-replay-versus-live fields are mandatory localized cross-representation ' ...
    'diagnostics; the fixed per-K table applies only to common-arithmetic parity'];
result.live_acceptance_classification = [ ...
    'the actual IBAMR solve must provide complete finite diagnostics and a fresh original ' ...
    'relative residual no larger than 1e-10'];
result.pass = result.candidate_original_residual_pass && ...
    result.candidate_physical_reporting_complete && ...
    result.oracle_physical_reporting_complete && ...
    result.iteration_history_reporting_complete;
end

function result = compare_histories(left, right)
left = left(:);
right = right(:);
result = struct('count_exact', numel(left) == numel(right), ...
                'left_count', int32(numel(left)), 'right_count', int32(numel(right)), ...
                'max_absolute_error', inf, 'denominator', 1, 'e_inf', inf);
if result.count_exact && ~isempty(left)
    left = left/max(abs(left(1)), 1e-30);
    right = right/max(abs(right(1)), 1e-30);
    comparison = compare_vectors(left, right);
    result.max_absolute_error = comparison.max_absolute_error;
    result.denominator = comparison.denominator;
    result.e_inf = comparison.e_inf;
end
end

function result = read_krylov_trace(filename)
contents = fileread(filename);
header = regexp(contents, '^ibamr-cav-krylov-trace-v1 ([0-9]+)$', ...
                'tokens', 'once', 'lineanchors');
if isempty(header)
    error('Unsupported Krylov trace schema in %s.', filename);
end
count = str2double(header{1});
tokens = regexp(contents, '^([0-9]+) ([0-9]+) ([^ ]+) ([^\n]+)$', ...
                'tokens', 'lineanchors');
if numel(tokens) ~= count
    error('Krylov trace record count mismatch in %s.', filename);
end
residuals = zeros(count, 1);
last_stem = '';
for k = 1:count
    if str2double(tokens{k}{1}) ~= k-1 || str2double(tokens{k}{2}) ~= k-1
        error('Noncontiguous Krylov trace in %s.', filename);
    end
    residuals(k) = str2double(tokens{k}{3});
    last_stem = strtrim(tokens{k}{4});
end
result = struct('count', count, 'residuals', residuals, 'last_stem', last_stem);
end

function result = read_named_scalars(filename)
contents = fileread(filename);
names = {'relative_residual', 'residual_momentum_l2', 'residual_divergence_l2', ...
         'ib_matrix_free_action_l2', 'ib_matrix_action_l2', 'ib_residual_l2', ...
         'ib_relative_residual'};
result = struct();
for k = 1:numel(names)
    name = names{k};
    token = regexp(contents, ['^' name ' ([^\n]+)$'], ...
                   'tokens', 'once', 'lineanchors');
    if isempty(token)
        error('Missing scalar %s in %s.', name, filename);
    end
    result.(name) = str2double(strtrim(token{1}));
end
end

function result = replay_policy(candidate_A, candidate_rhs, candidate_initial, candidate_actual, ...
                                candidate_patches, reference_patches, oracle_A, oracle_rhs, ...
                                pressure, target, directions)
direction_results = struct([]);
forward_state = [];
forward_fresh = [];
reverse_state = [];
reverse_fresh = [];
for k = 1:numel(directions)
    reverse = strcmp(directions{k}, 'reverse');
    current_direction = compare_sweeps(candidate_A, candidate_rhs, candidate_initial, ...
        candidate_patches, reference_patches, oracle_A, oracle_rhs, pressure, target, reverse);
    if isempty(direction_results)
        direction_results = current_direction;
    else
        direction_results(k) = current_direction;
    end
    if reverse
        reverse_state = direction_results(k).candidate_operator_final_state;
        reverse_fresh = direction_results(k).candidate_operator_final_residual;
    else
        forward_state = direction_results(k).candidate_operator_final_state;
        forward_fresh = direction_results(k).candidate_operator_final_residual;
    end
end

actual_live_control = struct( ...
    'correction', compare_vectors_gauged(forward_state, candidate_actual, pressure), ...
    'fresh_residual', compare_vectors(forward_fresh, candidate_rhs-candidate_A*candidate_actual));
actual_live_control.within_fixed_target = actual_live_control.correction.e_inf <= target && ...
                                          actual_live_control.fresh_residual.e_inf <= target;
actual_live_control.classification = ...
    'reported cross-runtime backend/arithmetic diagnostic; not a common-arithmetic gate';

if isempty(reverse_state)
    [reverse_state, reverse_fresh] = apply_sweep( ...
        candidate_A, candidate_rhs, candidate_initial, candidate_patches(end:-1:1), pressure);
end
reverse_sensitivity = struct( ...
    'correction', compare_vectors_gauged(forward_state, reverse_state, pressure), ...
    'fresh_residual', compare_vectors(forward_fresh, reverse_fresh));
reverse_sensitivity.detected = reverse_sensitivity.correction.e_inf > target || ...
                               reverse_sensitivity.fresh_residual.e_inf > target;

result = struct();
result.direction_results = strip_final_states(direction_results);
result.actual_live_forward_diagnostic_not_gate = actual_live_control;
result.reverse_order_sensitivity = reverse_sensitivity;
result.forward_state = forward_state;
result.pass = all([direction_results.pass]) && reverse_sensitivity.detected;
end

function results = strip_final_states(results)
results = rmfield(results, { ...
    'candidate_operator_final_state', 'candidate_operator_final_residual'});
end

function result = compare_sweeps(candidate_A, candidate_rhs, candidate_initial, ...
                                 candidate_patches, reference_patches, oracle_A, oracle_rhs, ...
                                 pressure, target, reverse)
if reverse
    direction = 'reverse';
    candidate_patches = candidate_patches(end:-1:1);
    reference_patches = reference_patches(end:-1:1);
else
    direction = 'forward';
end
[candidate_w, candidate_r, ~, candidate_trace] = apply_sweep( ...
    candidate_A, candidate_rhs, candidate_initial, candidate_patches, pressure);
[live_reference_w, live_reference_r, ~, live_reference_trace] = apply_sweep( ...
    candidate_A, candidate_rhs, candidate_initial, reference_patches, pressure);
[oracle_candidate_w, oracle_candidate_r, ~, oracle_candidate_trace] = apply_sweep( ...
    oracle_A, oracle_rhs, zeros(size(oracle_rhs)), candidate_patches, pressure);
[oracle_reference_w, oracle_reference_r, ~, oracle_reference_trace] = apply_sweep( ...
    oracle_A, oracle_rhs, zeros(size(oracle_rhs)), reference_patches, pressure);

result = struct();
result.direction = direction;
result.candidate_operator_every_patch = compare_sweep_traces( ...
    candidate_trace, live_reference_trace, target);
result.oracle_operator_every_patch = compare_sweep_traces( ...
    oracle_candidate_trace, oracle_reference_trace, target);
result.candidate_operator_final = paired_result( ...
    candidate_w, candidate_r, live_reference_w, live_reference_r, pressure, target);
result.oracle_operator_final = paired_result( ...
    oracle_candidate_w, oracle_candidate_r, oracle_reference_w, oracle_reference_r, pressure, target);
result.candidate_operator_final_state = candidate_w;
result.candidate_operator_final_residual = candidate_r;
result.pass = result.candidate_operator_every_patch.pass && ...
              result.oracle_operator_every_patch.pass && ...
              result.candidate_operator_final.pass && result.oracle_operator_final.pass;
end

function result = compare_sweep_traces(left, right, target)
if numel(left) ~= numel(right)
    error('Sweep traces have different patch counts.');
end
metrics = zeros(numel(left), 5);
for k = 1:numel(left)
    metrics(k, 1) = compare_matrices(left(k).A, right(k).A).e_inf;
    metrics(k, 2) = compare_vectors(left(k).rhs, right(k).rhs).e_inf;
    metrics(k, 3) = compare_vectors(left(k).delta, right(k).delta).e_inf;
    metrics(k, 4) = compare_vectors(left(k).state, right(k).state).e_inf;
    metrics(k, 5) = compare_vectors(left(k).residual, right(k).residual).e_inf;
end
first = find(any(metrics > target, 2), 1);
if isempty(first)
    first = -1;
else
    first = first-1;
end
result = struct( ...
    'patch_count', int32(numel(left)), ...
    'max_local_matrix_e_inf', max(metrics(:, 1)), ...
    'max_local_rhs_e_inf', max(metrics(:, 2)), ...
    'max_local_correction_e_inf', max(metrics(:, 3)), ...
    'max_global_state_e_inf', max(metrics(:, 4)), ...
    'max_residual_state_e_inf', max(metrics(:, 5)), ...
    'first_mismatch_traversal_step', int32(first), ...
    'pass', first < 0);
end

function value = read_candidate_state(directory, mapping, pressure, filename, equation_sign)
value = map_vector(read_vector_market(fullfile(directory, filename)), ...
                   mapping, pressure, equation_sign);
end

function result = read_native_diagnostics(directory, patch_count)
contents = fileread(fullfile(directory, 'candidate.stdout'));
result = struct( ...
    'sweep_count', read_output_number(contents, 'fgmres_local_sweep_count'), ...
    'solve_count', read_output_number(contents, 'fgmres_local_solve_count'), ...
    'max_backward_error', read_output_number(contents, 'fgmres_local_backward_error_max'), ...
    'max_incremental_fresh_error', read_output_number(contents, 'fgmres_incremental_fresh_error_max'), ...
    'max_local_rhs_error', read_output_number(contents, 'fgmres_local_rhs_error_max'), ...
    'diagnostics_valid', read_output_bool(contents, 'fgmres_local_diagnostics_valid'), ...
    'selected_trace_valid', read_output_bool(contents, 'cav_local_trace_selection_valid'));
result.expected_solve_count = result.sweep_count*patch_count;
result.pass = result.diagnostics_valid && result.selected_trace_valid && ...
              result.solve_count == result.expected_solve_count;
end

function result = validate_configuration(relaxed_dir, strict_dir, oracle_dir, ...
                                         stiffness, candidate_map, oracle_map)
relaxed_manifest = read_live_manifest(fullfile( ...
    relaxed_dir, 'cav_live_export_contract_manifest.txt'));
strict_manifest = read_live_manifest(fullfile( ...
    strict_dir, 'cav_live_export_contract_manifest.txt'));
oracle_metadata = jsondecode(fileread(fullfile(oracle_dir, 'oracle_metadata.json')));
relaxed_stdout = fileread(fullfile(relaxed_dir, 'candidate.stdout'));
strict_stdout = fileread(fullfile(strict_dir, 'candidate.stdout'));

result = struct();
result.dimension_2d = candidate_map.dimension == 2 && oracle_map.dimension == 2 && ...
                      oracle_metadata.dimension == 2;
result.global_dof_count = int32(candidate_map.count);
result.global_dof_count_exact = candidate_map.count == 192 && oracle_map.count == 192;
result.live_hierarchy_is_4_to_8 = oracle_metadata.base_grid_cells_per_axis == 4 && ...
    oracle_metadata.refinement_ratio == 2 && oracle_metadata.finest_grid_cells_per_axis == 8;
result.single_rank = str2double(relaxed_manifest.mpi_ranks) == 1 && ...
                     str2double(strict_manifest.mpi_ranks) == 1;
result.physical_stiffness_exact = ...
    read_output_number(relaxed_stdout, 'physical_stiffness') == stiffness && ...
    read_output_number(strict_stdout, 'physical_stiffness') == stiffness && ...
    oracle_metadata.physical_stiffness_K == stiffness;
result.relaxed_policy_exact = strcmp(relaxed_manifest.closure_policy, 'RELAXED');
result.strict_policy_exact = strcmp(strict_manifest.closure_policy, 'STRICT');
result.seed_and_order_exact = ...
    strcmp(relaxed_manifest.patch_seed_type, 'PRESSURE_CELL') && ...
    strcmp(strict_manifest.patch_seed_type, 'PRESSURE_CELL') && ...
    strcmp(relaxed_manifest.seed_stride, '1') && strcmp(strict_manifest.seed_stride, '1') && ...
    strcmp(relaxed_manifest.traversal_order, 'I_J') && strcmp(strict_manifest.traversal_order, 'I_J');
result.composition_and_backend_exact = ...
    strcmp(relaxed_manifest.composition, 'multiplicative') && ...
    strcmp(strict_manifest.composition, 'multiplicative') && ...
    strcmp(relaxed_manifest.local_solver_backend, 'blas-lapack-lu') && ...
    strcmp(strict_manifest.local_solver_backend, 'blas-lapack-lu');
result.sign_and_gauge_exact = ...
    strcmp(relaxed_manifest.pressure_equation, 'minus-div') && ...
    strcmp(strict_manifest.pressure_equation, 'minus-div') && ...
    str2double(relaxed_manifest.pressure_equation_row_multiplier_to_oracle) == -1 && ...
    str2double(strict_manifest.pressure_equation_row_multiplier_to_oracle) == -1 && ...
    strcmp(relaxed_manifest.pressure_gauge, 'zero-mean-correction') && ...
    strcmp(strict_manifest.pressure_gauge, 'zero-mean-correction') && ...
    strcmp(oracle_metadata.pressure_equation, 'div') && ...
    strcmp(oracle_metadata.pressure_gauge, 'zero-mean correction');
result.oracle_contract_exact = ...
    strcmp(oracle_metadata.oracle_sha, '5b77344db6746269f8c77695c99e9043907ba74b') && ...
    strcmp(oracle_metadata.patch_strategy, 'targeted_ib') && ...
    strcmp(oracle_metadata.composition, 'multiplicative') && ...
    strcmp(oracle_metadata.delta_function, 'IB_4') && oracle_metadata.patch_count == 64;
fields = fieldnames(result);
logical_fields = fields(~ismember(fields, {'global_dof_count'}));
result.pass = all(cellfun(@(name) logical(result.(name)), logical_fields));
end

function manifest = read_live_manifest(filename)
contents = fileread(filename);
if isempty(regexp(contents, '^ibamr-cav-live-export-v1$', 'once', 'lineanchors'))
    error('Unsupported live manifest schema in %s.', filename);
end
names = {'mpi_ranks', 'pressure_equation', ...
    'pressure_equation_row_multiplier_to_oracle', 'pressure_gauge', ...
    'patch_seed_type', 'closure_policy', 'seed_stride', 'traversal_order', ...
    'composition', 'local_solver_backend'};
manifest = struct();
for k = 1:numel(names)
    name = names{k};
    token = regexp(contents, ['^' name ' ([^\n]+)$'], ...
                   'tokens', 'once', 'lineanchors');
    if isempty(token)
        error('Missing live manifest field %s in %s.', name, filename);
    end
    manifest.(name) = strtrim(token{1});
end
end

function value = read_output_number(contents, name)
token = regexp(contents, ['^' name ' = ([^\n]+)$'], 'tokens', 'once', 'lineanchors');
if isempty(token)
    error('Missing native diagnostic %s.', name);
end
value = str2double(strtrim(token{1}));
end

function value = read_output_bool(contents, name)
token = regexp(contents, ['^' name ' = (true|false)$'], 'tokens', 'once', 'lineanchors');
if isempty(token)
    error('Missing native diagnostic %s.', name);
end
value = strcmp(token{1}, 'true');
end

function patches = construct_reference_patches(E_h, dof_map, pressure_seeds, policy)
velocity = strcmp(dof_map.kind, 'V');
pressure = strcmp(dof_map.kind, 'P');
N = round(sqrt(sum(pressure)));
if N^2 ~= sum(pressure)
    error('Reference patch construction requires a square 2D pressure grid.');
end

dof_by_key = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:dof_map.count
    dof_by_key(dof_key(dof_map, k)) = k;
end
adjacency = cell(dof_map.count, 1);
for row = find(velocity(:))'
    [~, columns, values] = find(E_h(row, :));
    row_max = max([0; abs(values(:))]);
    numerical_zero = max(numel(values)*eps*row_max, 1e-14*row_max);
    for entry = 1:numel(columns)
        column = columns(entry);
        if abs(values(entry)) <= numerical_zero
            continue;
        end
        if ~velocity(column)
            error('Reference E_h graph contains a pressure coupling.');
        end
        adjacency{row}(end+1) = column;
        adjacency{column}(end+1) = row;
    end
end
for k = 1:numel(adjacency)
    adjacency{k} = unique(adjacency{k});
end

patches = cell(numel(pressure_seeds), 1);
for seed_ordinal = 1:numel(pressure_seeds)
    seed = pressure_seeds(seed_ordinal);
    i = dof_map.i(seed);
    j = dof_map.j(seed);
    standard = standard_velocity_closure(dof_by_key, i, j, N);
    expanded = standard;
    for velocity_dof = standard(:)'
        expanded = union(expanded, adjacency{velocity_dof});
    end
    if isequal(expanded(:), standard(:))
        patches{seed_ordinal} = sort([standard(:); seed]);
        continue;
    end

    involved = seed;
    for velocity_dof = expanded(:)'
        axis = dof_map.axis(velocity_dof);
        vi = dof_map.i(velocity_dof);
        vj = dof_map.j(velocity_dof);
        if axis == 0
            involved(end+1, 1) = lookup_dof(dof_by_key, 'P', -1, mod(vi-1, N), vj); %#ok<AGROW>
            involved(end+1, 1) = lookup_dof(dof_by_key, 'P', -1, vi, vj); %#ok<AGROW>
        else
            involved(end+1, 1) = lookup_dof(dof_by_key, 'P', -1, vi, mod(vj-1, N)); %#ok<AGROW>
            involved(end+1, 1) = lookup_dof(dof_by_key, 'P', -1, vi, vj); %#ok<AGROW>
        end
    end
    involved = unique(involved);
    retained = [];
    for cell_dof = involved(:)'
        closure = standard_velocity_closure( ...
            dof_by_key, dof_map.i(cell_dof), dof_map.j(cell_dof), N);
        if strcmp(policy, 'RELAXED') || all(ismember(closure, expanded))
            retained(end+1, 1) = cell_dof; %#ok<AGROW>
        end
    end
    patch = retained;
    for cell_dof = retained(:)'
        patch = union(patch, standard_velocity_closure( ...
            dof_by_key, dof_map.i(cell_dof), dof_map.j(cell_dof), N));
    end
    if strcmp(policy, 'RELAXED')
        patch = union(patch, expanded);
    end
    patches{seed_ordinal} = sort(patch);
end
end

function dofs = standard_velocity_closure(dof_by_key, i, j, N)
dofs = [ ...
    lookup_dof(dof_by_key, 'V', 0, i, j); ...
    lookup_dof(dof_by_key, 'V', 0, mod(i+1, N), j); ...
    lookup_dof(dof_by_key, 'V', 1, i, j); ...
    lookup_dof(dof_by_key, 'V', 1, i, mod(j+1, N))];
dofs = sort(dofs);
end

function dof = lookup_dof(dof_by_key, kind, axis, i, j)
key = sprintf('%s:%d:%d:%d', kind, axis, i, j);
if ~isKey(dof_by_key, key)
    error('Missing reference DOF %s.', key);
end
dof = dof_by_key(key);
end

function label = stiffness_label(stiffness)
if stiffness == 1
    label = 'K1';
elseif stiffness == 1e2
    label = 'K1e2';
elseif stiffness == 1e4
    label = 'K1e4';
else
    error('Unsupported stiffness %.17g.', stiffness);
end
end

function target = fixed_target(stiffness)
if stiffness == 1
    target = 1e-13;
elseif stiffness == 1e2
    target = 1e-12;
elseif stiffness == 1e4
    target = 2e-11;
else
    error('Unsupported stiffness %.17g.', stiffness);
end
end

function pass = matrix_within(comparison, target)
pass = comparison.structure_exact && comparison.e_inf <= target;
end

function pass = matrix_exact(comparison)
pass = comparison.structure_exact && comparison.e_inf == 0;
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

function [w, fresh_residual, first_patch, trace] = apply_sweep(A, b, w, patches, pressure)
r = b-A*w;
first_patch = struct();
trace = repmat(struct('A', [], 'rhs', [], 'delta', [], 'state', [], 'residual', []), ...
               numel(patches), 1);
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
    trace(patch) = struct('A', local_A, 'rhs', local_rhs, 'delta', delta, ...
                          'state', w, 'residual', r);
end
w(pressure) = w(pressure)-mean(w(pressure));
fresh_residual = b-A*w;
end

function result = paired_result(left_w, left_r, right_w, right_r, pressure, target)
result = struct('correction', compare_vectors_gauged(left_w, right_w, pressure), ...
                'fresh_residual', compare_vectors(left_r, right_r));
result.pass = result.correction.e_inf <= target && result.fresh_residual.e_inf <= target;
end

function result = check_candidate_local_trace(candidate_dir, A, patch0_zero_based, target)
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
result.pass = result.matrix.e_inf <= target && ...
              result.restricted_rhs.e_inf <= target && ...
              result.direct_solve.e_inf <= target && result.backward_error <= target;
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
