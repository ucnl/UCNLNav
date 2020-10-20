% Solves the AOA navigation problem: searches for a horizontal angle of
% arrival by base points - points with known locations and measured times
% of arrival. Newton-Gauss method is used.

 function [a_rad, it_cnt] = Nav_aoa_ng_2d(bases, max_iterations, ...
     precision_threshold, el, velocity)
 
base_points_m = Nav_gcs_to_lcs(bases, el);

[a_rad, it_cnt] = Nav_tdoa_aoa_ng_2d_solve(base_points_m, max_iterations, ...
     precision_threshold, velocity);
 
end