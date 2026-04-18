function run_ex0_gmres_convergence_matlab_main()
  reference_root = getenv_default('CAV_REFERENCE_ROOT', '/Users/boyceg/code/implicit_ib_coupling_aware_vanka');
  summary_path = getenv_default('CAV_MATLAB_SUMMARY_PATH', '/tmp/matlab_summary.tsv');
  residual_dir = getenv_default('CAV_MATLAB_RESIDUAL_DIR', '/tmp/matlab_residual_histories');
  depths_csv = getenv_default('CAV_DEPTHS', '2,3,4');
  k_values_csv = getenv_default('CAV_K_VALUES', '1e2,1e4,1e6,1e8');
  coarse_n = str2double(getenv_default('CAV_COARSE_N', '4'));
  nu1 = str2double(getenv_default('CAV_NU1', '1'));
  nu2 = str2double(getenv_default('CAV_NU2', '1'));

  if isnan(coarse_n) || coarse_n < 2
    error('CAV_COARSE_N must be a positive integer >= 2');
  end
  if isnan(nu1) || nu1 < 1 || floor(nu1) ~= nu1
    error('CAV_NU1 must be a positive integer');
  end
  if isnan(nu2) || nu2 < 1 || floor(nu2) ~= nu2
    error('CAV_NU2 must be a positive integer');
  end

  addpath(reference_root);

  depths = parse_int_csv(depths_csv);
  k_tokens = parse_token_csv(k_values_csv);

  ensure_parent_dir(summary_path);
  ensure_dir(residual_dir);

  fid = fopen(summary_path, 'w');
  if fid < 0
    error('Unable to open summary path: %s', summary_path);
  end
  fprintf(fid, ['case_id\tdepth\tfinest_n\tcoarse_n\tK\tmatlab_converged\t' ...
                'matlab_iterations\tmatlab_relres\tmatlab_rhs_norm\tmatlab_dt\t' ...
                'matlab_rho\tmatlab_mu\tmatlab_lag_ds\tmatlab_num_curve_points\t' ...
                'matlab_nu1\tmatlab_nu2\tmatlab_history_file\n']);

  for depth_idx = 1:numel(depths)
    depth = depths(depth_idx);
    finest_n = coarse_n * (2^(depth - 1));

    for k_idx = 1:numel(k_tokens)
      k_token = strtrim(k_tokens{k_idx});
      if isempty(k_token)
        continue;
      end
      k_value = str2double(k_token);
      if isnan(k_value)
        error('Unable to parse K token: %s', k_token);
      end

      result = run_single_case(finest_n, depth, k_value, nu1, nu2);
      case_id = sprintf('depth%d_N%d_K%s', depth, finest_n, k_token);
      history_file = fullfile(residual_dir, [case_id, '.tsv']);
      write_history(history_file, result.resvec);

      fprintf(fid, ['%s\t%d\t%d\t%d\t%s\t%d\t%d\t%.17e\t%.17e\t%.17e\t' ...
                    '%.17e\t%.17e\t%.17e\t%d\t%d\t%d\t%s\n'], ...
              case_id, ...
              depth, ...
              finest_n, ...
              coarse_n, ...
              k_token, ...
              result.converged, ...
              result.iterations, ...
              result.relres, ...
              result.rhs_norm, ...
              result.dt, ...
              result.rho, ...
              result.mu, ...
              result.ds, ...
              result.num_curve_points, ...
              nu1, ...
              nu2, ...
              history_file);
      if exist('fflush', 'builtin') || exist('fflush', 'file')
        fflush(fid);
      end
    end
  end

  fclose(fid);
  fprintf(1, 'matlab_summary=%s\n', summary_path);
end

function result = run_single_case(N, depth, K, nu1, nu2)
  Lx = 1.0;
  Ly = 1.0;
  dx = Lx / N;
  dy = Ly / N;
  mfac = 1.0;
  r_cyl = 0.25;
  alpha = 0.23;
  beta = r_cyl^2 / alpha;
  rho = 1.0;
  mu = 1.0e-2;
  dt = 0.5 * dx;

  ds = mfac * dx / r_cyl;
  approx = round(2 * pi / ds);
  ds = 2 * pi / approx;
  l = length(0:ds:(2 * pi - ds));
  X_1 = Lx * (0.5 * ones(1, l) + alpha * cos(0:ds:(2 * pi - ds)));
  X_2 = Ly * (0.5 * ones(1, l) + beta * sin(0:ds:(2 * pi - ds)));
  X = [X_1; X_2]';

  zeros_u = zeros(N, N);
  zeros_v = zeros(N, N);
  Fn = apply_elastic_laplacian(X, K, ds);
  [fx_rhs, fy_rhs] = spreadIB4(Fn, X, zeros_u, zeros_v, dx, dy, ds);
  rhs = [fx_rhs(:); fy_rhs(:); zeros(N * N, 1)];

  L = sparse(build_saddle_point_matrix(X, N, N, dx, dy, ds, K, dt, rho, mu));

  box_size = 4;
  overlap = 2;
  relx_param = 1.0;
  M = @(x) v_cycle(L, x, zeros(size(x)), X, N, dx, ds, K, dt, rho, mu, ...
                   depth, nu1, nu2, box_size, overlap, relx_param, 1);

  restart = [];
  maxit = 150;
  tol = 1e-8;
  [~, flag, relres, ~, resvec] = gmres(L, rhs, restart, tol, maxit, M);

  result = struct();
  result.converged = int32(flag == 0);
  result.iterations = int32(max(0, length(resvec) - 1));
  result.relres = relres;
  result.resvec = resvec;
  result.rhs_norm = norm(rhs);
  result.dt = dt;
  result.rho = rho;
  result.mu = mu;
  result.ds = ds;
  result.num_curve_points = int32(approx);
end

function write_history(path, resvec)
  fid = fopen(path, 'w');
  if fid < 0
    error('Unable to open residual history file: %s', path);
  end
  for i = 1:length(resvec)
    fprintf(fid, '%d\t%.17e\n', i - 1, resvec(i));
  end
  fclose(fid);
end

function vals = parse_int_csv(raw)
  tokens = parse_token_csv(raw);
  vals = zeros(1, numel(tokens));
  for i = 1:numel(tokens)
    v = str2double(tokens{i});
    if isnan(v) || v < 1 || floor(v) ~= v
      error('Invalid integer token in CSV: %s', tokens{i});
    end
    vals(i) = int32(v);
  end
end

function tokens = parse_token_csv(raw)
  if isempty(raw)
    tokens = {};
    return;
  end
  parts = strsplit(raw, ',');
  tokens = {};
  for i = 1:numel(parts)
    t = strtrim(parts{i});
    if ~isempty(t)
      tokens{end + 1} = t; %#ok<AGROW>
    end
  end
end

function out = getenv_default(name, fallback)
  value = getenv(name);
  if isempty(value)
    out = fallback;
  else
    out = value;
  end
end

function ensure_parent_dir(path)
  [parent, ~, ~] = fileparts(path);
  if isempty(parent)
    return;
  end
  ensure_dir(parent);
end

function ensure_dir(path)
  if exist(path, 'dir') ~= 7
    [ok, msg] = mkdir(path);
    if ~ok
      error('Unable to create directory %s: %s', path, msg);
    end
  end
end

run_ex0_gmres_convergence_matlab_main();
