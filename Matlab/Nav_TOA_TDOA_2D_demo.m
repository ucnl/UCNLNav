% Nav_TOA_TDOA 2D_demo.m
close all;
clear all;
clc;


%% parameters
n_base_points = 4;
area_width_m = 1000;
max_depth_m = 100;
propagation_velocity = 1500;

max_iterations = 2000;
precision_threshold = 1E-9;
simplex_size = 1;

controur_levels = 32;

range_measurements_error = 0.01; % 0.01 means 1% of corresponding slant range

%% actual target location
r_ = rand * area_width_m / 2;
az_ = rand * 2 * pi;
actual_target_x = r_ * cos(az_);
actual_target_y = r_ * sin(az_);
actual_target_z = rand * max_depth_m;


%% base points
figure
hold on
grid on

base_points = zeros(n_base_points, 4);

for n = 1:n_base_points
   r_ = area_width_m / 2 - rand * area_width_m / 4;
   az_ = (n - 1) * 2 * pi / n_base_points;
   base_x = r_ * cos(az_);
   base_y = r_ * sin(az_);
   base_z = rand * max_depth_m;
   dist_to_target = Nav_dist_3d(base_x, base_y, base_z, actual_target_x, actual_target_y, actual_target_z);
   base_points(n, :) = [ base_x base_y base_z dist_to_target ];   
end

N =1:n_base_points;
plot3(actual_target_x, actual_target_y, actual_target_z,...
    'p',...
    'MarkerFaceColor', 'red',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

plot3(base_points(N, 1), base_points(N, 2), base_points(N, 3),...
    'o',...
    'MarkerFaceColor', 'green',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

for n = 1:n_base_points
   line([ actual_target_x, base_points(n, 1) ], [ actual_target_y, base_points(n, 2) ], [ actual_target_z, base_points(n, 3) ]); 
end

view(45, 15);
legend('target', 'base points');
title('Placement of base stations and the target');
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
zlabel('Z coordinate, m');


% adding range measurement errors
base_points(N, 4) = base_points(N, 4) + base_points(N, 4) *...
    (rand * range_measurements_error - range_measurements_error / 2);


% error surface tiles
tile_size_m = 10;
n_tiles = area_width_m / tile_size_m;


%% TOA solution
error_surface_toa = zeros(n_tiles, n_tiles);
for t_x = 1:n_tiles
   for t_y = 1:n_tiles      
      error_surface_toa(t_x, t_y) = Nav_eps_toa3d(base_points,...
          t_x * tile_size_m - area_width_m / 2,...
          t_y * tile_size_m - area_width_m / 2,...
          actual_target_z);
   end
end

figure
surf_a = [1:n_tiles] * tile_size_m - area_width_m / 2;
surf(surf_a, surf_a, error_surface_toa);
title('TOA solution: Residual function');
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
view(45, 15);

figure
hold on
contourf(surf_a, surf_a, error_surface_toa, controur_levels);
plot(actual_target_x, actual_target_y,...
    'p',...
    'MarkerFaceColor', 'red',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

plot(base_points(N, 1), base_points(N, 2),...
    'o',...
    'MarkerFaceColor', 'green',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

[ x_prev, y_prev ] = Nav_toa_circles_1d_solve(base_points, actual_target_z, pi / 180, 10, 0.1); 
[ x_best, y_best, rerr, it_cnt ] = Nav_toa_nlm_2d_solve(base_points, x_prev, y_prev, actual_target_z,...
    max_iterations, precision_threshold, simplex_size);

plot(x_best, y_best,...
    'd',...
    'MarkerFaceColor', 'yellow',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

title(sprintf('TOA Solution: Residual function. Target location estimated with E_{radial} = %.3f m in %d iterations', rerr, it_cnt));
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
legend('Residual function value', 'Actual target location', 'Base points', 'Estimated target location');


figure
hold on
grid on

dx = actual_target_x - x_best;
dy = actual_target_y - y_best;

plot(0, 0,...
    'p',...
    'MarkerFaceColor', 'red',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

plot(dx, dy,...
    'd',...
    'MarkerFaceColor', 'yellow',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(-dx * 2, -dy * 2, '.w');
plot(dx * 2, dy * 2, '.w');

d_delta = Nav_dist_3d(actual_target_x, actual_target_y, actual_target_z, x_best, y_best, actual_target_z);
title(sprintf('TOA Solution: Actual vs. Estimated location, distance: %.3f m', d_delta));
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
legend('Actual target location', 'Estimated target location');



%% TDOA solution

% since TDOA works with time difference of arriaval,
% we must recalculate base point's distances to times
base_points(N,4) = base_points(N,4) / propagation_velocity;
base_lines = Nav_build_base_lines(base_points, propagation_velocity);
error_surface_tdoa = zeros(n_tiles, n_tiles);
for t_x = 1:n_tiles
   for t_y = 1:n_tiles      
      error_surface_tdoa(t_x, t_y) = Nav_eps_tdoa3d(base_lines,...
          t_x * tile_size_m - area_width_m / 2,...
          t_y * tile_size_m - area_width_m / 2,...
          actual_target_z);
   end
end

figure
surf(surf_a, surf_a, error_surface_tdoa);
title('TDOA Solution: Residual function');
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
view(45, 15);

figure
hold on
contourf(surf_a, surf_a, error_surface_tdoa, controur_levels);
plot(actual_target_x, actual_target_y,...
    'p',...
    'MarkerFaceColor', 'red',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

plot(base_points(N, 1), base_points(N, 2),...
    'o',...
    'MarkerFaceColor', 'green',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

centroid = Nav_centroid(base_points);
[ x_best, y_best, rerr, it_cnt ] = Nav_tdoa_nlm_2d_solve(base_lines, centroid(1), centroid(2), actual_target_z,...
    max_iterations, precision_threshold, simplex_size);

plot(x_best, y_best,...
    'd',...
    'MarkerFaceColor', 'yellow',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

title(sprintf('TDOA Solution: Residual function. Target location estimated with E_{radial} = %.3f m in %d iterations', rerr, it_cnt));
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
legend('Residual function value', 'Actual target location', 'Base points', 'Estimated target location');


figure
hold on
grid on

dx = actual_target_x - x_best;
dy = actual_target_y - y_best;

plot(0, 0,...
    'p',...
    'MarkerFaceColor', 'red',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

plot(dx, dy,...
    'd',...
    'MarkerFaceColor', 'yellow',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(-dx * 2, -dy * 2, '.w');
plot(dx * 2, dy * 2, '.w');

d_delta = Nav_dist_3d(actual_target_x, actual_target_y, actual_target_z, x_best, y_best, actual_target_z);
title(sprintf('TDOA Solution: Actual vs. Estimated location, distance: %.3f m', d_delta));
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
legend('Actual target location', 'Estimated target location');


