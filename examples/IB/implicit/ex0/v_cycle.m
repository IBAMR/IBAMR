function [w, diagnostics] = v_cycle(L, b, w, X, N, h, ds, K, dt, rho, mu, depth, nu1, nu2, box_size, overlap, alpha, current_level)
% Parity override of MATLAB reference v_cycle() that additionally captures
% and exports per-level pre/post smoother solution snapshots.
% This override also routes subdomain extraction by closure policy:
% CAV_CLOSURE_POLICY=RELAXED -> extract_coupled_dofs()
% CAV_CLOSURE_POLICY=STRICT  -> extract_coupled_dofs_2()

    closure_policy = upper(strtrim(getenv('CAV_CLOSURE_POLICY')));
    if isempty(closure_policy)
        closure_policy = 'RELAXED';
    end
    if ~ismember(closure_policy, {'RELAXED', 'STRICT'})
        error('CAV_CLOSURE_POLICY must be RELAXED or STRICT');
    end

    if nargout > 1
        diagnostics = struct('captured', false, ...
                             'coarse_rhs', [], ...
                             'coarse_correction', [], ...
                             'pre_smooth_input', {cell(depth, 1)}, ...
                             'pre_smooth_output', {cell(depth, 1)}, ...
                             'post_smooth_input', {cell(depth, 1)}, ...
                             'post_smooth_output', {cell(depth, 1)});
    end

    if current_level == depth
        if nargout > 1
            diagnostics.captured = true;
            diagnostics.coarse_rhs = b;
        end
        w = L \ b;
        if nargout > 1
            diagnostics.coarse_correction = w;
        end
        return;
    end

    N_fine = N / (2^(current_level - 1));
    N_coarse = N_fine / 2;
    ib_level = depth - current_level;

    if nargout > 1
        current_level_pre_input = w;
    end
    for i = 1:nu1
        w = coupling_based_smoother_with_policy(L, b, w, N_fine, alpha, closure_policy);
    end
    if nargout > 1
        current_level_pre_output = w;
    end

    r = b - L * w;
    ru = reshape(r(1:N_fine^2), N_fine, N_fine);
    rv = reshape(r(N_fine^2 + 1:2 * N_fine^2), N_fine, N_fine);
    rp = reshape(r(2 * N_fine^2 + 1:3 * N_fine^2), N_fine, N_fine);
    [ru_coarse, rv_coarse] = apply_restriction_rt0_mac(ru, rv);
    rp_coarse = apply_pressure_restriction(rp);
    r_coarse = [ru_coarse(:); rv_coarse(:); rp_coarse(:)];

    L_coarse = build_saddle_point_matrix_level(X, N, h, ds, K, dt, rho, mu, current_level + 1);

    if nargout > 1
        [e_coarse, diagnostics] = v_cycle(L_coarse, r_coarse, zeros(size(r_coarse)), ...
            X, N, h, ds, K, dt, rho, mu, depth, nu1, nu2, box_size, overlap, alpha, current_level + 1);
        diagnostics.pre_smooth_input{ib_level + 1} = current_level_pre_input;
        diagnostics.pre_smooth_output{ib_level + 1} = current_level_pre_output;
    else
        e_coarse = v_cycle(L_coarse, r_coarse, zeros(size(r_coarse)), ...
            X, N, h, ds, K, dt, rho, mu, depth, nu1, nu2, box_size, overlap, alpha, current_level + 1);
    end

    [eu, ev] = apply_prolongation_rt0_mac( ...
        reshape(e_coarse(1:N_coarse^2), N_coarse, N_coarse), ...
        reshape(e_coarse(N_coarse^2 + 1:2 * N_coarse^2), N_coarse, N_coarse));
    ep = apply_pressure_prolongation( ...
        reshape(e_coarse(2 * N_coarse^2 + 1:3 * N_coarse^2), N_coarse, N_coarse));
    e = [eu(:); ev(:); ep(:)];

    w = w + e;

    if nargout > 1
        current_level_post_input = w;
    end
    for i = 1:nu2
        w = coupling_based_smoother_with_policy(L, b, w, N_fine, alpha, closure_policy);
    end
    if nargout > 1
        diagnostics.post_smooth_input{ib_level + 1} = current_level_post_input;
        diagnostics.post_smooth_output{ib_level + 1} = w;
        if current_level == 1
            export_smoother_snapshots(diagnostics, depth);
        end
    end
end

function w = coupling_based_smoother_with_policy(L, b, w, N, alpha, closure_policy)
% Policy-aware Vanka smoother used by parity-export v_cycle().
% RELAXED uses extract_coupled_dofs(); STRICT uses extract_coupled_dofs_2().

    n_u = N^2;
    stride = 1;
    for i = 1:stride:n_u
        r = b - L * w;
        if strcmp(closure_policy, 'STRICT')
            coupled_dofs = extract_coupled_dofs_2(L, i, N);
        else
            coupled_dofs = extract_coupled_dofs(L, i, N);
        end
        L_local = L(coupled_dofs, coupled_dofs);
        r_local = r(coupled_dofs);
        w_local = L_local \ r_local;
        w(coupled_dofs) = w(coupled_dofs) + alpha .* w_local;
    end

    % Maintain zero-mean pressure for periodic BC.
    w(2 * n_u + 1:end) = w(2 * n_u + 1:end) - mean(w(2 * n_u + 1:end));
end

function export_smoother_snapshots(diagnostics, depth)
    output_dir = getenv('CAV_MATLAB_PARITY_DIR');
    if isempty(output_dir)
        return;
    end

    for ln = 0:(depth - 1)
        idx = ln + 1;
        if idx <= numel(diagnostics.pre_smooth_input)
            v = diagnostics.pre_smooth_input{idx};
            if ~isempty(v)
                write_vector_matrix_market(fullfile(output_dir, sprintf('preconditioned_apply_pre_smooth_input_level%d.mtx', ln)), v);
            end
        end
        if idx <= numel(diagnostics.pre_smooth_output)
            v = diagnostics.pre_smooth_output{idx};
            if ~isempty(v)
                write_vector_matrix_market(fullfile(output_dir, sprintf('preconditioned_apply_pre_smooth_output_level%d.mtx', ln)), v);
                write_vector_matrix_market(fullfile(output_dir, sprintf('preconditioned_apply_pre_smooth_level%d.mtx', ln)), v);
            end
        end
        if idx <= numel(diagnostics.post_smooth_input)
            v = diagnostics.post_smooth_input{idx};
            if ~isempty(v)
                write_vector_matrix_market(fullfile(output_dir, sprintf('preconditioned_apply_post_smooth_input_level%d.mtx', ln)), v);
            end
        end
        if idx <= numel(diagnostics.post_smooth_output)
            v = diagnostics.post_smooth_output{idx};
            if ~isempty(v)
                write_vector_matrix_market(fullfile(output_dir, sprintf('preconditioned_apply_post_smooth_output_level%d.mtx', ln)), v);
                write_vector_matrix_market(fullfile(output_dir, sprintf('preconditioned_apply_post_smooth_level%d.mtx', ln)), v);
            end
        end
    end
end

function write_vector_matrix_market(path, vec)
    vec = full(vec(:));
    nnz_count = nnz(vec);
    fid = fopen(path, 'w');
    if fid < 0
        error('Unable to open MatrixMarket output path: %s', path);
    end
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real general\n');
    fprintf(fid, '%d %d %d\n', length(vec), 1, nnz_count);
    for i = 1:length(vec)
        value = vec(i);
        if value ~= 0.0
            fprintf(fid, '%d 1 %.17e\n', i, value);
        end
    end
    fclose(fid);
end
