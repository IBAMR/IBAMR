%% generate_tuna_geometry.m
% Generate CARANGIFORM (tuna-like) body geometry for IBAMR
%
% Creates a fusiform (streamlined) body characteristic of carangiform
% swimmers like tuna, mackerel, and sharks.
%
% Output: eel2d.vertex file (yes, still called eel2d for compatibility)

clear all; close all; clc;

fprintf('============================================\n');
fprintf('Generating CARANGIFORM (Tuna) Geometry\n');
fprintf('============================================\n');

%% Body Parameters
L = 1.0;                % Body length (m)
max_width = 0.12;       % Maximum width (12% of length - more robust!)
N_body = 128;           % Number of points along centerline

%% Tuna Body Profile
% Carangiform swimmers have fusiform (football-shaped) bodies
% Width profile: w(x) = w_max * 16 * X^2 * (1-X)^2
% This is a symmetric profile with maximum at X ≈ 0.5

x_body = linspace(0, L, N_body)';
X = x_body / L;  % Normalized position

% Fusiform body width profile (maximum at mid-body)
% Using a beta distribution-like function: 16*X^2*(1-X)^2
% This gives a smooth, streamlined shape
width = max_width * 16 * (X.^2) .* ((1-X).^2) / max(16 * (X.^2) .* ((1-X).^2));

% Add a small caudal peduncle (narrow tail base)
% Tuna have a characteristic narrow region before the tail
peduncle_start = 0.75;
peduncle_factor = 0.6;  % 60% of normal width
peduncle_mask = 1.0 - (1.0 - peduncle_factor) * ...
                exp(-50 * (X - peduncle_start).^2) .* (X > peduncle_start);
width = width .* peduncle_mask;

% Add small caudal fin (tail enlargement)
tail_start = 0.90;
tail_factor = 1.3;  % 30% wider
tail_mask = 1.0 + (tail_factor - 1.0) * ...
            exp(-100 * (X - tail_start - 0.08).^2) .* (X > tail_start);
width = width .* tail_mask;

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
filename = 'eel2d.vertex';  % Keep name for compatibility with input file
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
figure('Position', [100 100 1400 600]);

subplot(2,3,[1 4]);
plot(x_all, y_all, 'b-', 'LineWidth', 2);
hold on;
plot(x_all(1:10:end), y_all(1:10:end), 'ro', 'MarkerSize', 6);
axis equal; grid on;
xlabel('x (m)', 'FontSize', 12);
ylabel('y (m)', 'FontSize', 12);
title('CARANGIFORM (Tuna) Body Geometry', 'FontSize', 14, 'FontWeight', 'bold');
legend('Body outline', 'Lagrangian points', 'Location', 'best');

subplot(2,3,2);
plot(x_body, width, 'b-', 'LineWidth', 2);
hold on;
plot([0 L], [max_width max_width], 'r--', 'LineWidth', 1);
% Mark peduncle and tail regions
xline(peduncle_start*L, 'g--', 'Peduncle', 'LineWidth', 1.5);
xline(tail_start*L, 'm--', 'Caudal fin', 'LineWidth', 1.5);
grid on;
xlabel('Position along body (m)', 'FontSize', 12);
ylabel('Body width (m)', 'FontSize', 12);
title('Width Profile', 'FontSize', 14, 'FontWeight', 'bold');
legend('Actual width', 'Max width', 'Location', 'best');
xlim([0 L]);

subplot(2,3,3);
% Show amplitude envelope for comparison
X_plot = linspace(0, 1, 100);
a0 = 0.01; a1 = -0.05; a2 = 0.14;
A_envelope = a0 + a1*X_plot + a2*X_plot.^2;
plot(X_plot, A_envelope, 'r-', 'LineWidth', 2);
hold on;
plot(X_plot, width/max(width)*0.1, 'b--', 'LineWidth', 1.5);
grid on;
xlabel('Normalized position (X)', 'FontSize', 12);
ylabel('Amplitude / Width (normalized)', 'FontSize', 12);
title('Amplitude Envelope vs Body Width', 'FontSize', 12, 'FontWeight', 'bold');
legend('Amplitude A(X)', 'Body width', 'Location', 'best');

subplot(2,3,5);
% Fineness ratio along body
fineness = zeros(size(x_body));
for i = 1:length(x_body)
    if i == 1
        local_length = x_body(2) - x_body(1);
    elseif i == length(x_body)
        local_length = x_body(end) - x_body(end-1);
    else
        local_length = x_body(i+1) - x_body(i-1);
    end
    fineness(i) = local_length / (width(i) + 1e-10);
end
plot(x_body, fineness, 'g-', 'LineWidth', 2);
grid on;
xlabel('Position along body (m)', 'FontSize', 12);
ylabel('Local fineness ratio', 'FontSize', 12);
title('Streamlining (Length/Width)', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0 20]);

subplot(2,3,6);
% Cross-sectional area
area = pi * (width/2).^2;  % Assuming circular cross-section
plot(x_body, area*1e4, 'c-', 'LineWidth', 2);
grid on;
xlabel('Position along body (m)', 'FontSize', 12);
ylabel('Cross-sectional area (cm²)', 'FontSize', 12);
title('Body Cross-Section Area', 'FontSize', 12, 'FontWeight', 'bold');

sgtitle('Carangiform Swimming Geometry - Fusiform Body with Caudal Peduncle', ...
        'FontSize', 16, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'tuna_geometry.png');
fprintf('Geometry plot saved to: tuna_geometry.png\n');

%% Summary Statistics
fprintf('\n=== Geometry Summary ===\n');
fprintf('Body length: %.3f m\n', L);
fprintf('Max width: %.4f m (%.1f%%)\n', max_width, max_width/L*100);
fprintf('Aspect ratio: %.1f:1\n', L/max_width);
fprintf('Number of Lagrangian points: %d\n', N_total);
fprintf('Point spacing: ~%.4f m\n', L/(N_body-1));
fprintf('\nCharacteristics:\n');
fprintf('  - Fusiform (streamlined) body shape\n');
fprintf('  - Maximum width at ~50%% body length\n');
fprintf('  - Narrow caudal peduncle at 75%% length\n');
fprintf('  - Enlarged caudal fin region\n');
fprintf('  - Optimized for high-speed cruising\n');
fprintf('  - Suitable for carangiform swimming\n');
fprintf('============================================\n');
