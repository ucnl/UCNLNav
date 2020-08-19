% Nav_TOA_TDOA_3D_demo.m

close all;
clear all;
clc;


%% parameters
n_base_points = 4;
area_width_m = 1000;
max_depth_m = -100;
propagation_velocity = 1500;

max_iterations = 2000;
precision_threshold = 1E-12;
precision_threshold_hjs = 1E-4;
simplex_size = 3;
hjs_step_size = 3;

contour_levels = 32;

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
   base_z = -1.5;%rand * max_depth_m;
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
      error_surface_toa(t_y, t_x) = Nav_eps_toa3d(base_points,...
          (t_x - 1) * tile_size_m - area_width_m / 2,...
          (t_y - 1) * tile_size_m - area_width_m / 2,...
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
contourf(surf_a, surf_a, error_surface_toa, contour_levels);
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

[ x_best, y_best, z_best, rerr, it_cnt ] = Nav_toa_nlm_3d_solve(base_points, centroid(1), centroid(2), centroid(3),...
    max_iterations, precision_threshold, simplex_size);

[ x_best_hjs, y_best_hjs, z_best_hjs, rerr_hjs, it_cnt_hjs ] = Nav_toa_hjs_3d_solve(base_points, centroid(1), centroid(2), centroid(3),...
    max_iterations, precision_threshold_hjs, hjs_step_size);

plot(x_best, y_best,...
    'd',...
    'MarkerFaceColor', 'yellow',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(x_best_hjs, y_best_hjs,...
    'd',...
    'MarkerFaceColor', 'green',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

title('TOA Solution: Residual function.');
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
legend('Residual function value',...
       'Actual target location',...
       'Base points',...
       sprintf('Estimated target location (NLM), E_{radial} = %.3f m in %d iterations', rerr, it_cnt),...
       sprintf('Estimated target location (HJS), E_{radial} = %.3f m in %d iterations', rerr_hjs, it_cnt_hjs));

figure
hold on
grid on

dx_nlm = actual_target_x - x_best;
dy_nlm = actual_target_y - y_best;

dx_hjs = actual_target_x - x_best_hjs;
dy_hjs = actual_target_y - y_best_hjs;

dx = max(dx_nlm, dx_hjs);
dy = max(dy_nlm, dy_hjs);

plot(0, 0,...
    'p',...
    'MarkerFaceColor', 'red',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

plot(dx_nlm, dy_nlm,...
    'd',...
    'MarkerFaceColor', 'yellow',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(dx_hjs, dy_hjs,...
    'd',...
    'MarkerFaceColor', 'green',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(-dx * 2, -dy * 2, '.w');
plot(dx * 2, dy * 2, '.w');

d_delta_nlm = Nav_dist_3d(actual_target_x, actual_target_y, actual_target_z, x_best, y_best, actual_target_z);
d_delta_hjs = Nav_dist_3d(actual_target_x, actual_target_y, actual_target_z, x_best_hjs, y_best_hjs, actual_target_z);

title('TOA Solution: Actual vs. Estimated location');
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
legend('Actual target location',...
       sprintf('Estimated target location (NLM), \\sigma_{xy}: %.3f m, \\sigma_{z}: %.3f m', d_delta_nlm, actual_target_z - z_best),...
       sprintf('Estimated target location (HJS), \\sigma_{xy}: %.3f m, \\sigma_{z}: %.3f m', d_delta_hjs, actual_target_z - z_best_hjs));

   
%% 3D TOA - residual function along vertical axis

zmax = min([ z_best z_best_hjs max_depth_m]);
zmin = max([z_best z_best_hjs 0]);

zs = zmax:0.1:zmin;
toa3d_z_nlm = zeros(length(zs), 1);
toa3d_z_hjs = zeros(length(zs), 1);
for n = 1:length(zs)
    toa3d_z_nlm(n) = Nav_eps_toa3d(base_points, x_best, y_best, zs(n));
    toa3d_z_hjs(n) = Nav_eps_toa3d(base_points, x_best_hjs, y_best_hjs, zs(n));    
end

mn = min(toa3d_z_nlm);
mx = max(toa3d_z_nlm);
toa3d_z_nlm = (toa3d_z_nlm - mn) / (mx - mn);

mn = min(toa3d_z_hjs);
mx = max(toa3d_z_hjs);
toa3d_z_hjs = (toa3d_z_hjs - mn) / (mx - mn);

figure
hold on
grid on

plot(zs, toa3d_z_nlm, 'b');
plot(zs, toa3d_z_hjs, 'r');

title('Residual function along Z axis');
xlabel('Z coordinate');
ylabel('Normalized residual function value');

line([ actual_target_z actual_target_z ], [ 0 1 ], 'color', 'green');
line([ z_best z_best ], [ 0 1 ], 'color', 'blue');
line([ z_best_hjs z_best_hjs ], [ 0 1 ], 'color', 'red');
line([ max_depth_m max_depth_m ], [ 0 1 ], 'color', 'black');

zmaxh = zmax / 2;

if (max_depth_m < zmaxh)
    text(max_depth_m, 1, '\leftarrow Bottom depth');
else
    text(max_depth_m, 1, 'Bottom depth \rightarrow ','HorizontalAlignment','right')
end

if (z_best < zmaxh)
    text(z_best, 0.75, '\leftarrow Estimated Z (NLM)');
else
    text(z_best, 0.75, 'Estimated Z (NLM) \rightarrow ','HorizontalAlignment','right')
end

if (z_best_hjs < zmaxh)
    text(z_best_hjs, 0.5, '\leftarrow Estimated Z (HJS)');
else
    text(z_best_hjs, 0.5, 'Estimated Z (HJS) \rightarrow ','HorizontalAlignment','right')
end

if (actual_target_z < zmaxh)
    text(actual_target_z, 0.25, '\leftarrow Actual target Z');
else
    text(actual_target_z, 0.25, 'Actual target Z \rightarrow ','HorizontalAlignment','right')
end
    

%% TDOA solution

% since TDOA works with time difference of arriaval,
% we must recalculate base point's distances to times
base_points(N,4) = base_points(N,4) / propagation_velocity;
base_lines = Nav_build_base_lines(base_points, propagation_velocity);
error_surface_tdoa = zeros(n_tiles, n_tiles);
for t_x = 1:n_tiles
   for t_y = 1:n_tiles      
      error_surface_tdoa(t_y, t_x) = Nav_eps_tdoa3d(base_lines,...
          (t_x - 1) * tile_size_m - area_width_m / 2,...
          (t_y - 1) * tile_size_m - area_width_m / 2,...
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
contourf(surf_a, surf_a, error_surface_tdoa, contour_levels);
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

[ x_best, y_best, z_best, rerr, it_cnt ] = Nav_tdoa_nlm_3d_solve(base_lines, centroid(1), centroid(2), centroid(3),...
    max_iterations, precision_threshold, simplex_size);

[ x_best_hjs, y_best_hjs, z_best_hjs, rerr_hjs, it_cnt_hjs ] = Nav_tdoa_hjs_3d_solve(base_lines, centroid(1), centroid(2), centroid(3),...
    max_iterations, precision_threshold_hjs, hjs_step_size);


plot(x_best, y_best,...
    'd',...
    'MarkerFaceColor', 'yellow',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(x_best_hjs, y_best_hjs,...
    'd',...
    'MarkerFaceColor', 'green',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

title('TDOA Solution: Residual function.');
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
legend('Residual function value',...
       'Actual target location',...
       'Base points',...
       sprintf('Estimated target location (NLM), E_{radial} = %.3f m in %d iterations', rerr, it_cnt),...
       sprintf('Estimated target location (HJS), E_{radial} = %.3f m in %d iterations', rerr_hjs, it_cnt_hjs));


figure
hold on
grid on

dx_nlm = actual_target_x - x_best;
dy_nlm = actual_target_y - y_best;

dx_hjs = actual_target_x - x_best_hjs;
dy_hjs = actual_target_y - y_best_hjs;

dx = max(dx_nlm, dx_hjs);
dy = max(dy_nlm, dy_hjs);

plot(0, 0,...
    'p',...
    'MarkerFaceColor', 'red',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

plot(dx_nlm, dy_nlm,...
    'd',...
    'MarkerFaceColor', 'yellow',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(dx_hjs, dy_hjs,...
    'd',...
    'MarkerFaceColor', 'green',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(-dx * 2, -dy * 2, '.w');
plot(dx * 2, dy * 2, '.w');

d_delta_nlm = Nav_dist_3d(actual_target_x, actual_target_y, actual_target_z, x_best, y_best, actual_target_z);
d_delta_hjs = Nav_dist_3d(actual_target_x, actual_target_y, actual_target_z, x_best_hjs, y_best_hjs, actual_target_z);

title('TDOA Solution: Actual vs. Estimated location');
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
legend('Actual target location',...
       sprintf('Estimated target location (NLM), \\sigma_{xy}: %.3f m, \\sigma_{z}: %.3f m', d_delta_nlm, actual_target_z - z_best),...
       sprintf('Estimated target location (HJS), \\sigma_{xy}: %.3f m, \\sigma_{z}: %.3f m', d_delta_hjs, actual_target_z - z_best_hjs));


%% 3D TDOA - residual function along vertical axis

zmax = min([ z_best z_best_hjs max_depth_m]);
zmin = max([z_best z_best_hjs 0]);

zs = zmax:0.1:zmin;
tdoa3d_z_nlm = zeros(length(zs), 1);
tdoa3d_z_hjs = zeros(length(zs), 1);

for n = 1:length(zs)
    tdoa3d_z_nlm(n) = Nav_eps_tdoa3d(base_lines, x_best, y_best, zs(n));
    tdoa3d_z_hjs(n) = Nav_eps_tdoa3d(base_lines, x_best_hjs, y_best_hjs, zs(n));    
end

mn = min(tdoa3d_z_nlm);
mx = max(tdoa3d_z_nlm);
tdoa3d_z_nlm = (tdoa3d_z_nlm - mn) / (mx - mn);

mn = min(tdoa3d_z_hjs);
mx = max(tdoa3d_z_hjs);
tdoa3d_z_hjs = (tdoa3d_z_hjs - mn) / (mx - mn);

figure
hold on
grid on

plot(zs, tdoa3d_z_nlm, 'b');
plot(zs, tdoa3d_z_hjs, 'r');

title('Residual function along Z axis');
xlabel('Z coordinate');
ylabel('Normalized residual function value');

line([ actual_target_z actual_target_z ], [ 0 1 ], 'color', 'green');
line([ z_best z_best ], [ 0 1 ], 'color', 'blue');
line([ z_best_hjs z_best_hjs ], [ 0 1 ], 'color', 'red');
line([ max_depth_m max_depth_m ], [ 0 1 ], 'color', 'black');

zmaxh = zmax / 2;

if (max_depth_m < zmaxh)
    text(max_depth_m, 1, '\leftarrow Bottom depth');
else
    text(max_depth_m, 1, 'Bottom depth \rightarrow ','HorizontalAlignment','right')
end

if (z_best < zmaxh)
    text(z_best, 0.75, '\leftarrow Estimated Z (NLM)');
else
    text(z_best, 0.75, 'Estimated Z (NLM) \rightarrow ','HorizontalAlignment','right')
end

if (z_best_hjs < zmaxh)
    text(z_best_hjs, 0.5, '\leftarrow Estimated Z (HJS)');
else
    text(z_best_hjs, 0.5, 'Estimated Z (HJS) \rightarrow ','HorizontalAlignment','right')
end

if (actual_target_z < zmaxh)
    text(actual_target_z, 0.25, '\leftarrow Actual target Z');
else
    text(actual_target_z, 0.25, 'Actual target Z \rightarrow ','HorizontalAlignment','right')
end


%% 3D TDOA problem: trying to estimate at least x and y

figure
hold on
contourf(surf_a, surf_a, error_surface_tdoa, contour_levels);
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

[ x_best, y_best, rerr, it_cnt ] = Nav_tdoa_nlm_2d_solve(base_lines, centroid(1), centroid(2), centroid(3),...
    max_iterations, precision_threshold, simplex_size);

[ x_best_hjs, y_best_hjs, rerr_hjs, it_cnt_hjs ] = Nav_tdoa_hjs_2d_solve(base_lines, centroid(1), centroid(2), centroid(3),...
    max_iterations, precision_threshold_hjs, hjs_step_size);


plot(x_best, y_best,...
    'd',...
    'MarkerFaceColor', 'yellow',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(x_best_hjs, y_best_hjs,...
    'd',...
    'MarkerFaceColor', 'green',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

title('2D TDOA Solution on 3D TDOA Problem: Residual function.');
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
legend('Residual function value',...
       'Actual target location',...
       'Base points',...
       sprintf('Estimated target location (NLM), E_{radial} = %.3f m in %d iterations', rerr, it_cnt),...
       sprintf('Estimated target location (HJS), E_{radial} = %.3f m in %d iterations', rerr_hjs, it_cnt_hjs));


figure
hold on
grid on

dx_nlm = actual_target_x - x_best;
dy_nlm = actual_target_y - y_best;

dx_hjs = actual_target_x - x_best_hjs;
dy_hjs = actual_target_y - y_best_hjs;

dx = max(dx_nlm, dx_hjs);
dy = max(dy_nlm, dy_hjs);

plot(0, 0,...
    'p',...
    'MarkerFaceColor', 'red',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 15);

plot(dx_nlm, dy_nlm,...
    'd',...
    'MarkerFaceColor', 'yellow',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(dx_hjs, dy_hjs,...
    'd',...
    'MarkerFaceColor', 'green',...
    'MarkerEdgeColor', 'blue',...
    'MarkerSize', 7);

plot(-dx * 2, -dy * 2, '.w');
plot(dx * 2, dy * 2, '.w');

d_delta_nlm = Nav_dist_3d(actual_target_x, actual_target_y, actual_target_z, x_best, y_best, actual_target_z);
d_delta_hjs = Nav_dist_3d(actual_target_x, actual_target_y, actual_target_z, x_best_hjs, y_best_hjs, actual_target_z);

title('2D TDOA Solution on 3D TDOA Problem: Actual vs. Estimated location');
xlabel('X coordinate, m');
ylabel('Y coordinate, m');
legend('Actual target location',...
       sprintf('Estimated target location (NLM), \\sigma_{xy}: %.3f m', d_delta_nlm),...
       sprintf('Estimated target location (HJS), \\sigma_{xy}: %.3f m', d_delta_hjs));
