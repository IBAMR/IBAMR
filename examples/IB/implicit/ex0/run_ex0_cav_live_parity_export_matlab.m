function run_ex0_cav_live_parity_export_matlab_main()
  wrapper_dir = fileparts(mfilename('fullpath'));
  matlab_reference_repo = getenv_default('CAV_MATLAB_REFERENCE_REPO', '/Users/boyceg/code/implicit_ib_coupling_aware_vanka');
  reference_script = fullfile(matlab_reference_repo, 'run_cav_live_parity_export.m');
  override_v_cycle = fullfile(wrapper_dir, 'v_cycle.m');
  original_path = path;
  local_scratch_dir = tempname(tempdir);
  local_reference_copy = fullfile(local_scratch_dir, 'run_ex0_cav_live_parity_export_reference_copy.m');
  local_v_cycle_copy = fullfile(local_scratch_dir, 'v_cycle.m');

  if ~isfolder(matlab_reference_repo)
    error('CAV_MATLAB_REFERENCE_REPO does not exist: %s', matlab_reference_repo);
  end
  if ~isfile(reference_script)
    error('MATLAB reference export script not found: %s', reference_script);
  end
  if ~isfile(override_v_cycle)
    error('MATLAB parity override v_cycle.m not found: %s', override_v_cycle);
  end
  [ok, msg] = mkdir(local_scratch_dir);
  if ~ok
    error('Unable to create temporary scratch directory %s: %s', local_scratch_dir, msg);
  end
  cleanup = onCleanup(@() cleanup_wrapper(original_path, local_scratch_dir)); %#ok<NASGU>

  copyfile(reference_script, local_reference_copy, 'f');
  copyfile(override_v_cycle, local_v_cycle_copy, 'f');
  restoredefaultpath();
  addpath(local_scratch_dir, '-begin');
  addpath(matlab_reference_repo, '-end');
  run(local_reference_copy);
end

function cleanup_wrapper(original_path, local_scratch_dir)
  path(original_path);
  if exist(local_scratch_dir, 'dir') == 7
    [ok, msg] = rmdir(local_scratch_dir, 's');
    if ~ok
      warning('Failed to remove temporary scratch directory %s: %s', local_scratch_dir, msg);
    end
  end
end

function value = getenv_default(name, fallback)
  value = getenv(name);
  if isempty(value)
    value = fallback;
  end
end

run_ex0_cav_live_parity_export_matlab_main();
