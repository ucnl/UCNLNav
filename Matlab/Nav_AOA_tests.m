% Nav_AOA_tests.m
close all;
clear all;%variables;
clc;

%% more parameters
range_m = 1000;
propagation_velocity = 1500;
base_size = 1.5;

n_base_points = 4;
max_iterations = 1000;
precision_threshold = 1E-12;
contour_levels = 32;
time_measurement_error = 0.01; % 0.01 means 1% of corresponding slant range
uncorrelated_time_measurement_error = 0.0001;

base_points = zeros(n_base_points, 4);

for n = 1:n_base_points   
   az_ = (n - 1) * 2 * pi / (n_base_points);
   base_x = base_size * cos(az_);
   base_y = base_size * sin(az_);
   base_z = 0;
   base_points(n, :) = [ base_x base_y base_z 0 ];   
end

errors = [];
error_idx = 1;

for actual_angle = 0:359
   
   actual_angle_rad = actual_angle * pi / 180;
      
   target_x = range_m * cos(actual_angle_rad);
   target_y = range_m * sin(actual_angle_rad);
   target_z = 0 * 100;
          
   for n = 1:n_base_points
       base_points(n, 4) = Nav_dist_3d(base_points(n, 1), base_points(n, 2), base_points(n, 3),...
           target_x, target_y, target_z) / propagation_velocity;              
   end
   
   mean_time = mean(base_points(:,4));
   common_error = mean_time * rand * time_measurement_error;
   
   for n = 1:n_base_points
       base_points(n, 4) = base_points(n, 4) + common_error + rand * uncorrelated_time_measurement_error;
   end
      
   [a_rad, it_cnt] = Nav_tdoa_aoa_ng_2d_solve(base_points, max_iterations,...
     precision_threshold, propagation_velocity);
    
   errors(error_idx) = actual_angle_rad - a_rad;
   if (errors(error_idx) > pi)
      errors(error_idx) = 2 * pi - errors(error_idx); 
   end
   
   error_idx = error_idx + 1;
end


figure
plot(errors * 180 / pi);
title(sprintf('AOA estimation errors with random measurements noise (uncorrelated: %.5f sec, correlated: %.1f %% of slant range)', uncorrelated_time_measurement_error, time_measurement_error * 100));
xlabel('Actual angle of arrival, °');
ylabel('Angle of arrival estimation error, °');
