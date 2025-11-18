%% ========================================================================
%% NACA0012 ANGUILLIFORM (EEL-LIKE) MESH GENERATOR
%% Based on IBAMR eel2d example (eel2d_straightswimmer.m)
%% Reference: Flow transitions and mapping for undulating swimmers, Muhammad Saif Ullah Khalid
%% ========================================================================

clear all; clc;

fprintf('========================================================================\n');
fprintf('NACA0012 ANGUILLIFORM (EEL-LIKE) MESH GENERATOR\n');
fprintf('========================================================================\n\n');

%% ========================================================================
%% PARAMETERS
%% ========================================================================

% Background Mesh Properties (match eel2d for comparison)
Ny = 16*4*4*4;     % 1024 grid points in y
Ly = 4;            % Domain height (4L)
Nx = 32*4*4*4;     % 2048 grid points in x
Lx = 8;            % Domain length (8L)

% Calculate grid spacing
dy = Ly/Ny;        % ~0.00390625
dx = Lx/Nx;        % ~0.00390625

% Fish geometry
L = 1;             % Fish body length (non-dimensional)

% Anguilliform amplitude envelope coefficients
% A(x/L) = 0.0367 + 0.0323(x/L) + 0.0310(x/L)^2
% Gradual amplitude increase from head to tail (eel-like)
a0 = 0.0367;       % Constant term (head amplitude)
a1 = 0.0323;       % Linear coefficient (gradual increase)
a2 = 0.0310;       % Quadratic coefficient (tail acceleration)

% NACA0012 thickness parameters
thickness_scale = 0.2489;  % Scale to match eel2d point density
t_max = 0.12 * thickness_scale;  % Effective thickness ~4.2%

% Rotation (0 for horizontal swimming)
angleRotation = 0;
RotationMatrix = [cos(angleRotation) -sin(angleRotation);
                  sin(angleRotation)  cos(angleRotation)];

fprintf('Parameters:\n');
fprintf('  Domain: %.0fL × %.0fL\n', Lx, Ly);
fprintf('  Grid: %d × %d\n', Nx, Ny);
fprintf('  Spacing: dx = dy = %.10f\n', dx);
fprintf('  Amplitude: A(x/L) = %.4f + %.4f(x/L) + %.4f(x/L)^2\n', a0, a1, a2);
fprintf('  NACA thickness: %.1f%% (scaled)\n\n', t_max*100);

%% ========================================================================
%% MESH GENERATION
%% ========================================================================

fprintf('Generating Lagrangian mesh...\n');

BodyNx = ceil(L/dx);  % Number of cross-sections
Coord{BodyNx,2} = []; % Storage for coordinates
NumMatPoints = 0;     % Counter for total points

% Generate cross-sections along body
for i = 1:BodyNx

    x = (i-1)*dx;        % Streamwise position
    x_norm = x/L;        % Normalized position (0 to 1)

    % Centerline position: y(x,t=0) = A(x/L) * cos(2π*x/L)
    % For anguilliform, amplitude increases gradually along entire body
    y = (a0 + a1*x_norm + a2*x_norm^2) * cos(2*pi*x);

    % NACA0012 thickness distribution
    y_t = 5*t_max*(0.2969*sqrt(x_norm) - 0.1260*x_norm - ...
                   0.3516*x_norm^2 + 0.2843*x_norm^3 - 0.1015*x_norm^4);

    % Full thickness (both surfaces)
    height = 2*y_t*L;

    % Number of points across thickness
    NumPointsInHeight = max(1, ceil(height/dy));

    % Initialize arrays
    ycoord_up = []; ycoord_down = [];
    xcoord_up = []; xcoord_down = [];

    % Generate vertical line of points
    for j = 1:NumPointsInHeight

        xshifted = x - 0.5;  % Center at origin
        yshifted = y;        % On deformed centerline

        RotatedCoord = RotationMatrix*[xshifted; yshifted];

        xcoord_up(j)   = RotatedCoord(1) - (j-1)*dy*sin(angleRotation);
        xcoord_down(j) = RotatedCoord(1) + j*dy*sin(angleRotation);
        ycoord_up(j)   = RotatedCoord(2) + (j-1)*dy*cos(angleRotation);
        ycoord_down(j) = RotatedCoord(2) - j*dy*cos(angleRotation);
        NumMatPoints = NumMatPoints + 2;
    end

    Coord{i,1} = cat(2, xcoord_up, xcoord_down);
    Coord{i,2} = cat(2, ycoord_up, ycoord_down);
end

fprintf('✓ Mesh complete: %d IB points\n\n', NumMatPoints);

%% ========================================================================
%% VISUALIZATION 1: AMPLITUDE ENVELOPE
%% ========================================================================

fprintf('Generating amplitude envelope plot...\n');

x_plot = linspace(0, 1, 200);
A_plot = a0 + a1*x_plot + a2*x_plot.^2;

figure('Position', [100, 100, 800, 500], 'Color', 'w');
plot(x_plot, A_plot, 'b-', 'LineWidth', 3);
hold on;
plot([0 1], [a0 a0+a1+a2], 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

set(gca, 'FontSize', 14, 'FontName', 'Times');
xlabel('x/L', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('A', 'FontSize', 16, 'FontWeight', 'bold');
title('NACA0012 Anguilliform Amplitude Envelope', 'FontSize', 16, 'FontWeight', 'bold');

text(0.05, a0+0.005, sprintf('Head: A(0) = %.4f', a0), 'FontSize', 12);
text(0.85, a0+a1+a2+0.005, sprintf('Tail: A(1) = %.4f', a0+a1+a2), 'FontSize', 12);

grid on; xlim([0 1]); ylim([0 0.12]); box on;

print('amplitude_envelope_anguilliform', '-dpng', '-r300');
fprintf('✓ Saved: amplitude_envelope_anguilliform.png\n');

%% ========================================================================
%% VISUALIZATION 2: BACKBONE MOTION
%% ========================================================================

fprintf('Generating backbone motion plot...\n');

figure('Position', [100, 100, 1000, 500], 'Color', 'w');

x_backbone = linspace(0, 1, 100);
num_phases = 10;

for phase_idx = 1:num_phases
    t_phase = (phase_idx-1) / num_phases;
    A_backbone = a0 + a1*x_backbone + a2*x_backbone.^2;
    y_backbone = A_backbone .* cos(2*pi*(x_backbone - t_phase));

    color_val = 0.2 + 0.6*t_phase;
    plot(x_backbone, y_backbone, '-', 'LineWidth', 2, ...
         'Color', [color_val, 0.2, 1-color_val]);
    hold on;
end

set(gca, 'FontSize', 14, 'FontName', 'Times');
xlabel('x/L', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('y/L', 'FontSize', 16, 'FontWeight', 'bold');
title('NACA0012 Anguilliform Backbone Motion (10 phases)', ...
      'FontSize', 16, 'FontWeight', 'bold');

grid on; xlim([0 1]); ylim([-0.12 0.12]); box on;

text(0.5, 0.09, 'y(x,t) = A(x/L) × cos[2π(x/L - ft)]', ...
     'FontSize', 13, 'HorizontalAlignment', 'center', ...
     'BackgroundColor', 'w', 'EdgeColor', 'k');

print('backbone_motion_anguilliform', '-dpng', '-r300');
fprintf('✓ Saved: backbone_motion_anguilliform.png\n');

%% ========================================================================
%% VISUALIZATION 3: MESH
%% ========================================================================

fprintf('Generating mesh plot...\n');

figure('Position', [100, 100, 1000, 500], 'Color', 'w');

for i = 1:size(Coord,1)
    plot(Coord{i,1}(:), Coord{i,2}(:), '.', 'MarkerSize', 2, ...
         'Color', [0 0.45 0.74]);
    hold on;
end

set(gca, 'FontSize', 14, 'FontName', 'Times');
xlabel('x/L', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('y/L', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('NACA0012 Anguilliform Mesh (%d IB points)', NumMatPoints), ...
      'FontSize', 16, 'FontWeight', 'bold');

grid on; axis equal; box on;

text(-0.45, 0.05, sprintf('Grid: %d × %d\nSpacing: %.6f\nDomain: %.0fL × %.0fL', ...
     Nx, Ny, dx, Lx, Ly), 'FontSize', 11, ...
     'BackgroundColor', 'w', 'EdgeColor', 'k');

print('mesh_visualization_anguilliform', '-dpng', '-r300');
fprintf('✓ Saved: mesh_visualization_anguilliform.png\n');

%% ========================================================================
%% COMPARISON PLOT: ANGUILLIFORM VS CARANGIFORM ENVELOPES
%% ========================================================================

fprintf('Generating mode comparison plot...\n');

% Anguilliform envelope
A_anguilli = a0 + a1*x_plot + a2*x_plot.^2;

% Carangiform envelope (from naca_carangiform_clean.m)
a0_caran = 0.02;
a1_caran = -0.0825;
a2_caran = 0.1625;
A_carangi = a0_caran + a1_caran*x_plot + a2_caran*x_plot.^2;

figure('Position', [100, 100, 900, 600], 'Color', 'w');
plot(x_plot, A_anguilli, 'b-', 'LineWidth', 3, 'DisplayName', 'Anguilliform (Eel-like)');
hold on;
plot(x_plot, A_carangi, 'r-', 'LineWidth', 3, 'DisplayName', 'Carangiform (Fish-like)');

% Mark key points
plot(0, a0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(1, a0+a1+a2, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(0, a0_caran, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(1, a0_caran+a1_caran+a2_caran, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

set(gca, 'FontSize', 14, 'FontName', 'Times');
xlabel('x/L (Normalized Position)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Amplitude A(x/L)', 'FontSize', 16, 'FontWeight', 'bold');
title('Swimming Mode Comparison: Amplitude Envelopes', 'FontSize', 16, 'FontWeight', 'bold');

legend('Location', 'best', 'FontSize', 13);
grid on; xlim([0 1]); ylim([0 0.12]); box on;

% Add annotations
text(0.5, 0.105, 'Anguilliform: Gradual increase', 'FontSize', 11, ...
     'Color', 'b', 'HorizontalAlignment', 'center');
text(0.75, 0.085, 'Carangiform: Posterior concentration', 'FontSize', 11, ...
     'Color', 'r', 'HorizontalAlignment', 'center');

print('swimming_modes_comparison', '-dpng', '-r300');
fprintf('✓ Saved: swimming_modes_comparison.png\n');

%% ========================================================================
%% WRITE VERTEX FILE
%% ========================================================================

fprintf('\nWriting vertex file...\n');

filename = 'naca0012anguilliform.vertex';
fid = fopen(filename, 'wt');
fprintf(fid, '%d\n', NumMatPoints);

for i = 1:size(Coord,1)
    for j = 1:length(Coord{i,1})
        fprintf(fid, '%f\t%f\n', Coord{i,1}(j), Coord{i,2}(j));
    end
end

fclose(fid);
fprintf('✓ Saved: %s\n', filename);

%% ========================================================================
%% SUMMARY
%% ========================================================================

fprintf('\n========================================================================\n');
fprintf('COMPLETE!\n');
fprintf('========================================================================\n');
fprintf('Swimming Mode:      ANGUILLIFORM (Eel-like)\n');
fprintf('IB points:          %d\n', NumMatPoints);
fprintf('Cross-sections:     %d\n', BodyNx);
fprintf('Points/section:     %.1f (average)\n', NumMatPoints/BodyNx);
fprintf('Amplitude:          A(0) = %.4f (head), A(1) = %.4f (tail)\n', ...
        a0, a0+a1+a2);
fprintf('\nFiles generated:\n');
fprintf('  1. naca0012anguilliform.vertex (IBAMR input)\n');
fprintf('  2. amplitude_envelope_anguilliform.png\n');
fprintf('  3. backbone_motion_anguilliform.png\n');
fprintf('  4. mesh_visualization_anguilliform.png\n');
fprintf('  5. swimming_modes_comparison.png\n');
fprintf('========================================================================\n');
fprintf('\nTo use this mesh in IBAMR:\n');
fprintf('  1. Update input2d: structure_names = "naca0012anguilliform"\n');
fprintf('  2. Uncomment anguilliform equations in ConstraintIBKinematics section\n');
fprintf('  3. Comment out carangiform equations\n');
fprintf('========================================================================\n');
