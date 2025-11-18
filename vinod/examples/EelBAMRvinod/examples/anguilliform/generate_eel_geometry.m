%% generate_eel_geometry.m
% Generate ANGUILLIFORM (eel-like) body geometry for IBAMR
%
% Creates a slender, gradually tapering body characteristic of anguilliform
% swimmers like eels and lampreys.
%
% Output: eel2d.vertex file

clear all; close all; clc;

fprintf('============================================\n');
fprintf('Generating ANGUILLIFORM (Eel) Geometry\n');
fprintf('============================================\n');

%% Body Parameters
L = 1.0;                % Body length (m)
max_width = 0.04;       % Maximum width (4% of length - slender!)
N_body = 128;           % Number of points along centerline

%% Eel Body Profile
% Anguilliform swimmers have slender, gradually tapering bodies
% Width profile: w(x) = w_max * sin(πx/L)^0.7

x_body = linspace(0, L, N_body)';
X = x_body / L;  % Normalized position

% Slender body width profile (slightly thicker in middle)
width = max_width * (sin(pi * X).^0.7);

%% Generate upper and lower surfaces
x_upper = x_body;
y_upper = width / 2;

x_lower = flipud(x_body);
y_lower = -flipud(width / 2);

%% Combine into closed contour
x_all = [x_upper; x_lower];
y_all = [y_upper; y_lower];
N_total = length(x_all);

fprintf('Total points: %d\n', N_total);
fprintf('Max width: %.4f m (%.1f%% of length)\n', max_width, max_width/L*100);

%% Center the fish at origin
x_all = x_all - L/2;

%% Write vertex file for IBAMR
filename = 'eel2d.vertex';
fid = fopen(filename, 'w');

% Write header: number of vertices
fprintf(fid, '%d\n', N_total);

% Write vertex coordinates
for i = 1:N_total
    fprintf(fid, '%.12f %.12f\n', x_all(i), y_all(i));
end

fclose(fid);
fprintf('Geometry saved to: %s\n', filename);

%% Visualization
figure('Position', [100 100 1200 400]);

subplot(1,2,1);
plot(x_all, y_all, 'b-', 'LineWidth', 2);
hold on;
plot(x_all(1:10:end), y_all(1:10:end), 'ro', 'MarkerSize', 6);
axis equal; grid on;
xlabel('x (m)', 'FontSize', 12);
ylabel('y (m)', 'FontSize', 12);
title('ANGUILLIFORM (Eel) Body Geometry', 'FontSize', 14, 'FontWeight', 'bold');
legend('Body outline', 'Lagrangian points', 'Location', 'best');

subplot(1,2,2);
plot(x_body, width, 'b-', 'LineWidth', 2);
hold on;
plot([0 L], [max_width max_width], 'r--', 'LineWidth', 1);
grid on;
xlabel('Position along body (m)', 'FontSize', 12);
ylabel('Body width (m)', 'FontSize', 12);
title('Width Profile', 'FontSize', 14, 'FontWeight', 'bold');
legend('Actual width', 'Max width', 'Location', 'best');
xlim([0 L]);

sgtitle('Anguilliform Swimming Geometry - Slender Body', 'FontSize', 16, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'eel_geometry.png');
fprintf('Geometry plot saved to: eel_geometry.png\n');

%% Summary Statistics
fprintf('\n=== Geometry Summary ===\n');
fprintf('Body length: %.3f m\n', L);
fprintf('Max width: %.4f m (%.1f%%)\n', max_width, max_width/L*100);
fprintf('Aspect ratio: %.1f:1\n', L/max_width);
fprintf('Number of Lagrangian points: %d\n', N_total);
fprintf('Point spacing: ~%.4f m\n', L/(N_body-1));
fprintf('\nCharacteristics:\n');
fprintf('  - Slender, gradually tapering body\n');
fprintf('  - High aspect ratio (typical of eels)\n');
fprintf('  - Smooth width distribution\n');
fprintf('  - Suitable for anguilliform swimming\n');
fprintf('============================================\n');
