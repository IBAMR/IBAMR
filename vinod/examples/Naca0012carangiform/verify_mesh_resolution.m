%% ========================================================================
%% MESH RESOLUTION VERIFICATION SCRIPT
%% Checks consistency between MATLAB mesh generator and IBAMR input2d
%% ========================================================================

clear all; clc;

fprintf('========================================================================\n');
fprintf('NACA0012 MESH RESOLUTION VERIFICATION\n');
fprintf('========================================================================\n\n');

%% ========================================================================
%% SECTION 1: MATLAB Mesh Generator Parameters
%% ========================================================================

fprintf('1. MATLAB MESH GENERATOR PARAMETERS\n');
fprintf('   ----------------------------------------\n');

% Your current parameters
Ny_matlab = 16*4*4*4;     % 1024 grid points in y
Ly = 4;                   % Domain height (4L)
Nx_matlab = 32*4*4*4;     % 2048 grid points in x
Lx = 8;                   % Domain length (8L)
L = 1;                    % Fish body length

% Calculate grid spacing
dy_matlab = Ly/Ny_matlab;
dx_matlab = Lx/Nx_matlab;

fprintf('   Nx (x-direction):     %d\n', Nx_matlab);
fprintf('   Ny (y-direction):     %d\n', Ny_matlab);
fprintf('   Lx (domain length):   %.2f L\n', Lx);
fprintf('   Ly (domain height):   %.2f L\n', Ly);
fprintf('   dx (grid spacing):    %.10f\n', dx_matlab);
fprintf('   dy (grid spacing):    %.10f\n', dy_matlab);
fprintf('   Isotropic:            %s\n\n', ...
        iif(abs(dx_matlab-dy_matlab)<1e-10, 'YES ✓', 'NO ✗'));

%% ========================================================================
%% SECTION 2: IBAMR input2d Parameters
%% ========================================================================

fprintf('2. IBAMR INPUT2D CONFIGURATION\n');
fprintf('   ----------------------------------------\n');

% From input2d file
N_input2d = 64;           % Coarsest grid cells in x
MAX_LEVELS = 3;           % Maximum AMR levels (0,1,2,3)
REF_RATIO = 4;            % Refinement ratio

% Domain from input2d
x_lo = [-6.52, -1.52];
x_up = [1.48, 2.48];
Lx_input2d = x_up(1) - x_lo(1);
Ly_input2d = x_up(2) - x_lo(2);

% Calculate finest level resolution
Nx_finest = 2*N_input2d * REF_RATIO^MAX_LEVELS;
Ny_finest = N_input2d * REF_RATIO^MAX_LEVELS;
dx_finest = Lx_input2d / Nx_finest;
dy_finest = Ly_input2d / Ny_finest;

fprintf('   N (coarse x-cells):   %d\n', N_input2d);
fprintf('   MAX_LEVELS:           %d (levels 0-%d)\n', MAX_LEVELS, MAX_LEVELS);
fprintf('   REF_RATIO:            %d\n', REF_RATIO);
fprintf('   Lx (from x_lo/x_up):  %.2f\n', Lx_input2d);
fprintf('   Ly (from x_lo/x_up):  %.2f\n\n', Ly_input2d);

fprintf('   FINEST AMR LEVEL (Level %d):\n', MAX_LEVELS);
fprintf('   Nx (finest):          %d\n', Nx_finest);
fprintf('   Ny (finest):          %d\n', Ny_finest);
fprintf('   dx (finest):          %.10f\n', dx_finest);
fprintf('   dy (finest):          %.10f\n', dy_finest);
fprintf('   Isotropic:            %s\n\n', ...
        iif(abs(dx_finest-dy_finest)<1e-10, 'YES ✓', 'NO ✗'));

%% ========================================================================
%% SECTION 3: Comparison and Consistency Check
%% ========================================================================

fprintf('3. CONSISTENCY CHECK\n');
fprintf('   ----------------------------------------\n');

% Check domain sizes
domain_match = (abs(Lx - Lx_input2d) < 0.01) && (abs(Ly - Ly_input2d) < 0.01);
fprintf('   Domain size match:    %s\n', iif(domain_match, 'YES ✓', 'NO ✗'));

% Check resolution
Nx_match = (Nx_matlab == Nx_finest);
Ny_match = (Ny_matlab == Ny_finest);
resolution_match = Nx_match && Ny_match;

fprintf('   Nx match:             %s ', iif(Nx_match, 'YES ✓', 'NO ✗'));
if ~Nx_match
    fprintf('(MATLAB: %d, input2d: %d, ratio: %.2f)', ...
            Nx_matlab, Nx_finest, Nx_finest/Nx_matlab);
end
fprintf('\n');

fprintf('   Ny match:             %s ', iif(Ny_match, 'YES ✓', 'NO ✗'));
if ~Ny_match
    fprintf('(MATLAB: %d, input2d: %d, ratio: %.2f)', ...
            Ny_matlab, Ny_finest, Ny_finest/Ny_matlab);
end
fprintf('\n');

% Check grid spacing
dx_match = abs(dx_matlab - dx_finest) < 1e-8;
dy_match = abs(dy_matlab - dy_finest) < 1e-8;
spacing_match = dx_match && dy_match;

fprintf('   dx match:             %s ', iif(dx_match, 'YES ✓', 'NO ✗'));
if ~dx_match
    fprintf('(MATLAB: %.8f, input2d: %.8f, ratio: %.2f)', ...
            dx_matlab, dx_finest, dx_finest/dx_matlab);
end
fprintf('\n');

fprintf('   dy match:             %s ', iif(dy_match, 'YES ✓', 'NO ✗'));
if ~dy_match
    fprintf('(MATLAB: %.8f, input2d: %.8f, ratio: %.2f)', ...
            dy_matlab, dy_finest, dy_finest/dy_matlab);
end
fprintf('\n\n');

%% ========================================================================
%% SECTION 4: Lagrangian Mesh Statistics
%% ========================================================================

fprintf('4. LAGRANGIAN MESH RESOLUTION\n');
fprintf('   ----------------------------------------\n');

% Calculate for current MATLAB parameters
BodyNx_matlab = ceil(L/dx_matlab);
fprintf('   Using MATLAB spacing (dx=%.6f):\n', dx_matlab);
fprintf('   - Cross-sections:     %d\n', BodyNx_matlab);

% Estimate IB points (very rough - actual depends on thickness)
est_points_per_section = 11.5;  % From your vertex file
est_total_points_matlab = BodyNx_matlab * est_points_per_section;
fprintf('   - Est. IB points:     ~%.0f\n\n', est_total_points_matlab);

% Calculate for input2d finest level
BodyNx_finest = ceil(L/dx_finest);
est_total_points_finest = BodyNx_finest * est_points_per_section;
fprintf('   Using input2d finest spacing (dx=%.6f):\n', dx_finest);
fprintf('   - Cross-sections:     %d\n', BodyNx_finest);
fprintf('   - Est. IB points:     ~%.0f\n\n', est_total_points_finest);

% Check if vertex file exists
vertex_file = 'naca0012carangiform.vertex';
if exist(vertex_file, 'file')
    fid = fopen(vertex_file, 'r');
    actual_points = fscanf(fid, '%d', 1);
    fclose(fid);
    fprintf('   Actual vertex file:   %s\n', vertex_file);
    fprintf('   - Actual IB points:   %d\n\n', actual_points);
else
    fprintf('   Vertex file not found: %s\n\n', vertex_file);
end

%% ========================================================================
%% SECTION 5: Recommendations
%% ========================================================================

fprintf('5. RECOMMENDATIONS\n');
fprintf('   ----------------------------------------\n');

if resolution_match && spacing_match
    fprintf('   ✓ Configuration is CONSISTENT!\n');
    fprintf('   Your MATLAB mesh generator parameters match the finest\n');
    fprintf('   AMR level in input2d. No changes needed.\n\n');
else
    fprintf('   ⚠ RESOLUTION MISMATCH DETECTED!\n\n');

    fprintf('   The MATLAB mesh generator uses different spacing than\n');
    fprintf('   the finest AMR level in input2d.\n\n');

    fprintf('   OPTION A: Keep MATLAB parameters, adjust input2d\n');
    fprintf('   ------------------------------------------------\n');
    fprintf('   To match Nx=%d, Ny=%d at finest level:\n', Nx_matlab, Ny_matlab);

    % Calculate required MAX_LEVELS
    % Nx_matlab = 2*N*REF_RATIO^MAX_LEVELS
    % Solve for MAX_LEVELS
    if Nx_matlab == 2048 && Ny_matlab == 1024
        fprintf('   Change input2d to:\n');
        fprintf('     N = 128\n');
        fprintf('     MAX_LEVELS = 2  (levels 0,1,2)\n');
        fprintf('     REF_RATIO = 4\n');
        fprintf('   This gives finest: 128*4^2*2 = 2048 (x), 128*4^2 = 1024 (y)\n\n');
    else
        fprintf('   Custom calculation needed for your Nx,Ny values\n\n');
    end

    fprintf('   OPTION B: Keep input2d, adjust MATLAB\n');
    fprintf('   ------------------------------------------------\n');
    fprintf('   Change MATLAB mesh generator to:\n');
    fprintf('     Nx = %d\n', Nx_finest);
    fprintf('     Ny = %d\n', Ny_finest);
    fprintf('     dx = %.10f\n', dx_finest);
    fprintf('     dy = %.10f\n', dy_finest);
    fprintf('   Then regenerate .vertex file\n\n');

    fprintf('   OPTION C: Use intermediate resolution\n');
    fprintf('   ------------------------------------------------\n');
    fprintf('   Use MAX_LEVELS-1 as target resolution:\n');
    Nx_intermediate = 2*N_input2d * REF_RATIO^(MAX_LEVELS-1);
    Ny_intermediate = N_input2d * REF_RATIO^(MAX_LEVELS-1);
    dx_intermediate = Lx_input2d / Nx_intermediate;
    dy_intermediate = Ly_input2d / Ny_intermediate;
    fprintf('     Nx = %d\n', Nx_intermediate);
    fprintf('     Ny = %d\n', Ny_intermediate);
    fprintf('     dx = %.10f\n', dx_intermediate);
    fprintf('     dy = %.10f\n\n', dy_intermediate);
end

%% ========================================================================
%% SECTION 6: IB Resolution Quality Metrics
%% ========================================================================

fprintf('6. IB RESOLUTION QUALITY\n');
fprintf('   ----------------------------------------\n');

% Thickness parameters
thickness_scale = 0.2489;
t_max = 0.12 * thickness_scale;
max_thickness = 2 * 5 * t_max * L;  % Maximum body thickness

fprintf('   MATLAB configuration:\n');
points_across_matlab = max_thickness / dy_matlab;
fprintf('   - Max body thickness:     %.6f L\n', max_thickness);
fprintf('   - Points across body:     %.1f\n', points_across_matlab);
fprintf('   - Resolution quality:     %s\n', ...
        iif(points_across_matlab >= 8, 'EXCELLENT ✓', ...
        iif(points_across_matlab >= 4, 'GOOD ✓', 'POOR ✗')));

fprintf('\n   input2d finest configuration:\n');
points_across_finest = max_thickness / dy_finest;
fprintf('   - Points across body:     %.1f\n', points_across_finest);
fprintf('   - Resolution quality:     %s\n\n', ...
        iif(points_across_finest >= 8, 'EXCELLENT ✓', ...
        iif(points_across_finest >= 4, 'GOOD ✓', 'POOR ✗')));

% IB kernel support
fprintf('   IB_4 kernel requirements:\n');
fprintf('   - Needs ~4 grid points per IB wavelength\n');
fprintf('   - Lagrangian spacing should be ≈ Eulerian spacing\n');
fprintf('   - Your setup: ds_lag ≈ dx_matlab = %.6f  %s\n\n', ...
        dx_matlab, iif(abs(dx_matlab-dy_matlab)<1e-8, '✓', '✗'));

%% ========================================================================
%% SECTION 7: Computational Cost Estimate
%% ========================================================================

fprintf('7. COMPUTATIONAL COST ESTIMATE\n');
fprintf('   ----------------------------------------\n');

% Rough cost scaling: O(N_grid * N_lagrangian)
cost_matlab = Nx_matlab * Ny_matlab * est_total_points_matlab;
cost_finest = Nx_finest * Ny_finest * est_total_points_finest;
cost_ratio = cost_finest / cost_matlab;

fprintf('   MATLAB configuration:\n');
fprintf('   - Eulerian points:    %d × %d = %d\n', ...
        Nx_matlab, Ny_matlab, Nx_matlab*Ny_matlab);
fprintf('   - Lagrangian points:  ~%.0f\n', est_total_points_matlab);
fprintf('   - Relative cost:      1.0×\n\n');

fprintf('   input2d finest configuration:\n');
fprintf('   - Eulerian points:    %d × %d = %d\n', ...
        Nx_finest, Ny_finest, Nx_finest*Ny_finest);
fprintf('   - Lagrangian points:  ~%.0f\n', est_total_points_finest);
fprintf('   - Relative cost:      %.1f×\n\n', cost_ratio);

%% ========================================================================
%% SECTION 8: Summary
%% ========================================================================

fprintf('========================================================================\n');
fprintf('SUMMARY\n');
fprintf('========================================================================\n');

if resolution_match && spacing_match
    fprintf('✓ Configurations are consistent\n');
    fprintf('✓ Ready to run simulation\n');
else
    fprintf('⚠ Resolution mismatch detected\n');
    fprintf('→ Choose Option A, B, or C above to resolve\n');
end

fprintf('\nCurrent status:\n');
fprintf('  MATLAB mesh: %d × %d grid, dx=dy=%.6f\n', ...
        Nx_matlab, Ny_matlab, dx_matlab);
fprintf('  input2d finest: %d × %d grid, dx=dy=%.6f\n', ...
        Nx_finest, Ny_finest, dx_finest);
fprintf('  Lagrangian: ~%.0f IB points, %d cross-sections\n', ...
        est_total_points_matlab, BodyNx_matlab);
fprintf('========================================================================\n');

%% Helper function
function result = iif(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
