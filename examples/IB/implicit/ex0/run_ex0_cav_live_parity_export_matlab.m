function run_ex0_cav_live_parity_export_matlab_main()
  wrapper_dir = fileparts(mfilename('fullpath'));
  matlab_reference_repo = getenv_default('CAV_MATLAB_REFERENCE_REPO', '/Users/boyceg/code/implicit_ib_coupling_aware_vanka');
  reference_script = fullfile(matlab_reference_repo, 'run_cav_live_parity_export.m');
  local_scratch_dir = tempdir;
  local_reference_copy = fullfile(local_scratch_dir, 'run_ex0_cav_live_parity_export_reference_copy.m');
  local_v_cycle_copy = fullfile(local_scratch_dir, 'v_cycle.m');

  if ~isfolder(matlab_reference_repo)
    error('CAV_MATLAB_REFERENCE_REPO does not exist: %s', matlab_reference_repo);
  end
  if ~isfile(reference_script)
    error('MATLAB reference export script not found: %s', reference_script);
  end

  copyfile(reference_script, local_reference_copy, 'f');
  copyfile(fullfile(wrapper_dir, 'v_cycle.m'), local_v_cycle_copy, 'f');
  addpath(matlab_reference_repo);
  addpath(local_scratch_dir, '-begin');
  run(local_reference_copy);
end

function value = getenv_default(name, fallback)
  value = getenv(name);
  if isempty(value)
    value = fallback;
  end
end

run_ex0_cav_live_parity_export_matlab_main();
