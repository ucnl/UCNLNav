% Nav_tests.m
close all;
clear all;%variables;
clc;

VNC_DEF_IT_LIMIT = 2000;    % Default iterations limit for Vincenty algorithms
VNC_DEF_EPSILON = 1E-12;    % Default precision threshold for Vincenty algorithms
NLM_DEF_IT_LIMIT = 1200;    % Default iterations limit for Nelder-Mead optimization algorithm
NLM_DEF_PREC_THRLD = 1E-8;  % Default precision threshold for Nelder-Mead optimization algorithm

AOA_DEF_IT_LIMIT = 100;     % Default iterations limit for Newton-Gauss in AOA estimation algorithm
AOA_DEF_PREC_THRLD = 1E-12; % Default precision threshold for Newton-Gauss in AOA estimation algorithm

propagation_velocity = 1500.0; % in m/s
el = Nav_build_standard_ellipsoid('WGS84');



sp_lat_deg = 48.073566;
sp_lon_deg = -123.02983;
ep_lat_deg = 48.064458;
ep_lon_deg = -123.02983;

sp_lat_rad = degtorad(sp_lat_deg);
sp_lon_rad = degtorad(sp_lon_deg);
ep_lat_rad = degtorad(ep_lat_deg);
ep_lon_rad = degtorad(ep_lon_deg);

adist_m0 = Nav_haversine_inverse(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad, el.mjsa_m);
adist_m100 = Nav_haversine_inverse(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad, el.mjsa_m - 1000);

adist_m0 - adist_m100


%%
tic;

%%
fprintf('[ Testing ] Nav_BuildStandardEllipsoid... ');
assert_approx_eq(el.mjsa_m, 6378137.0, 1E-6);
assert_approx_eq(el.ifltn, 298.257223563, 1E-6);
assert_approx_eq(el.fltn, 1 / el.ifltn, 1E-6);
assert_approx_eq(el.mnsa_m, el.mjsa_m * (1 - el.fltn), 1E-6);
assert_approx_eq(el.ecnt, (el.mjsa_m^2 - el.mnsa_m^2) / el.mjsa_m^2, 1E-6);
fprintf('Ok\n');
    
%%
fprintf('[ Testing ] Nav_wrap...');
assert_approx_eq(Nav_wrap(10.0, 5.0), 5.0, 1E-6);
assert_approx_eq(Nav_wrap(-10.0, 5.0), -5.0, 1E-6);
assert_approx_eq(Nav_wrap(0.0, 5.0), 0.0, 1E-6);
assert_approx_eq(Nav_wrap(5.0, 10.0), 5.0, 1E-6);
assert_approx_eq(Nav_wrap(-5.0, 10.0), -5.0, 1E-6);
fprintf('Ok\n');

%%
fprintf('[ Testing ] Nav_lat_1deg_length...');
assert_approx_eq(Nav_lat_1deg_length(-3.141592654 / 2.0, el), Nav_lat_1deg_length(3.141592654 / 2.0, el), 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-3.141592654, el), 111314.502041079, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(3.141592654, el), 111314.502041079, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-3.141592654, el), Nav_lat_1deg_length(3.141592654, el), 1E-6);
assert_approx_eq(Nav_lat_1deg_length(0.000000000, el), 111314.502041079, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(0.000000000, el), 111314.502041079, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(0.174532925, el), 111314.727675276, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-0.174532925, el), 111314.727675276, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(0.349065850, el), 111315.377367309, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-0.349065850, el), 111315.377367309, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(0.523598776, el), 111316.372765512, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-0.523598776, el), 111316.372765512, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(0.698131701, el), 111317.593822430, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-0.698131701, el), 111317.593822430, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(0.872664626, el), 111318.893268579, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-0.872664626, el), 111318.893268579, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(1.047197551, el), 111320.114371577, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-1.047197551, el), 111320.114371577, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(1.221730476, el), 111321.109840379, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-1.221730476, el), 111321.109840379, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(1.396263402, el), 111321.759594497, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-1.396263402, el), 111321.759594497, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(1.570796327, el), 111321.985253213, 1E-6);
assert_approx_eq(Nav_lat_1deg_length(-1.570796327, el), 111321.985253213, 1E-6);

for n = -90:90
    n_rad = degtorad(n);
    assert_approx_eq(Nav_lat_1deg_length(n_rad, el), Nav_lat_1deg_length(-n_rad, el), 1E-3);
end
fprintf('Ok\n');

%%
fprintf('[ Testing ] Nav_lon_1deg_length...');
assert_approx_eq(Nav_lon_1deg_length(-3.141592654, el), 111319.490793274, 1E-5);
assert_approx_eq(Nav_lon_1deg_length(3.141592654, el), 111319.490793274, 1E-5);
assert_approx_eq(Nav_lon_1deg_length(-3.141592654, el), Nav_lon_1deg_length(3.141592654, el), 1E-5);
assert_approx_eq(Nav_lon_1deg_length(0.000000000, el), 111319.490793274, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(0.174532925, el), 109628.371666625, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-0.174532925, el), 109628.371666625, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(0.349065850, el), 104606.378238853, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-0.349065850, el), 104606.378238853, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(0.523598776, el), 96406.047016128, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-0.523598776, el), 96406.047016128, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(0.698131701, el), 85276.466841735, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-0.698131701, el), 85276.466841735, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(0.872664626, el), 71555.730303868, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-0.872664626, el), 71555.730303868, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(1.047197551, el), 55660.680811254, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-1.047197551, el), 55660.680811254, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(1.221730476, el), 38074.261548399, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-1.221730476, el), 38074.261548399, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(1.396263402, el), 19330.846811735, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-1.396263402, el), 19330.846811735, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(1.570796327, el), 0.000000000, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-1.570796327, el), 0.000000000, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(1.745329252, el), 19330.846811735, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-1.745329252, el), 19330.846811735, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(1.919862177, el), 38074.261548399, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-1.919862177, el), 38074.261548399, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(2.094395102, el), 55660.680811254, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-2.094395102, el), 55660.680811254, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(2.268928028, el), 71555.730303868, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-2.268928028, el), 71555.730303868, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(2.443460953, el), 85276.466841735, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-2.443460953, el), 85276.466841735, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(2.617993878, el), 96406.047016128, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-2.617993878, el), 96406.047016128, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(2.792526803, el), 104606.378238853, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-2.792526803, el), 104606.378238853, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(2.967059728, el), 109628.371666625, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-2.967059728, el), 109628.371666625, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(3.141592654, el), 111319.490793274, 1E-4);
assert_approx_eq(Nav_lon_1deg_length(-3.141592654, el), 111319.490793274, 1E-4);

for n = -180:180 
    n_rad = degtorad(n);
    assert_approx_eq(Nav_lon_1deg_length(n_rad, el), Nav_lon_1deg_length(-n_rad, el), 1E-5);
end
fprintf('Ok\n');

%%
fprintf('[ Testing ] Nav_geopoint_offset_by_deltas and Nav_get_deltas_by_geopoints...');
for lat_deg = -8:8
    lat_rad = degtorad(lat_deg * 10.0);
    for lon_deg = -17:17
       lon_rad = degtorad(lon_deg * 10.0);
       for delta_m = -1:3
           d_m = 10^delta_m;
           [ new_lat_rad, new_lon_rad ] = Nav_geopoint_offset_by_deltas(lat_rad, lon_rad, delta_m, delta_m, el);
           [ delta_lat_m, delta_lon_m ] = Nav_get_deltas_by_geopoints(lat_rad, lon_rad, new_lat_rad, new_lon_rad, el);
           assert_approx_eq(delta_lat_m, delta_m, 1E-4 * d_m);
           assert_approx_eq(delta_lon_m, delta_m, 1E-4 * d_m);
       end
    end      
end
fprintf('Ok\n');

%%
fprintf('[ Testing ] Nav_geopoint_offset_by_deltas_wgs84 and Nav_get_deltas_by_geopoints_wgs84...');
for lat_deg = -8:8
    lat_rad = degtorad(lat_deg * 10.0);
    for lon_deg = -17:17
       lon_rad = degtorad(lon_deg * 10.0);
       for delta_m = -1:3
           d_m = 10^delta_m;
           [ new_lat_rad, new_lon_rad ] = Nav_geopoint_offset_by_deltas_wgs84(lat_rad, lon_rad, delta_m, delta_m);
           [ delta_lat_m, delta_lon_m ] = Nav_get_deltas_by_geopoints_wgs84(lat_rad, lon_rad, new_lat_rad, new_lon_rad);
           assert_approx_eq(delta_lat_m, delta_m, 1E-4 * d_m);
           assert_approx_eq(delta_lon_m, delta_m, 1E-4 * d_m);
       end
    end      
end
fprintf('Ok\n');

%%
fprintf('[ Testing ] Nav_haversine_inverse, Nav_haversine_direct and Nav_haversine_initial_bearing...');
for lat_deg = -8:8
    sp_lat_rad = degtorad(lat_deg * 10.0);
    for lon_deg = -17:17  
        sp_lon_rad = degtorad(lon_deg * 10.0);
        for az_deg = 0:35  
            fwd_az_rad = Nav_wrap_2pi(az_deg * 10.0 * pi / 180);  
            for dist_m = -1:3     
                dst_m = 10.0^dist_m;    
                [ ep_lat_rad, ep_lon_rad ] = Nav_haversine_direct(sp_lat_rad, sp_lon_rad, dst_m, fwd_az_rad, el.mjsa_m);
                adist_m = Nav_haversine_inverse(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad, el.mjsa_m);
                a_az_rad = Nav_haversine_initial_bearing(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad);
                if abs(abs(a_az_rad - fwd_az_rad) - 2 * pi) < 10E-3 
                    assert_approx_eq(abs(a_az_rad - fwd_az_rad), 2 * pi, 10E-3);                               
                else 
                    assert_approx_eq(fwd_az_rad, a_az_rad, 10E-3);
                end
                assert_approx_eq(dst_m, adist_m, 10E-1);
            end
        end
    end
end
fprintf('Ok\n');

%%
fprintf('[ Testing ] Nav_vincenty_inverse and Nav_vincenty_direct...');   
for lat_deg = -9:9
    sp_lat_rad = Nav_wrap_2pi(degtorad(lat_deg * 10.0));
    for lon_deg = -18:18
        sp_lon_rad = Nav_wrap_2pi(degtorad(lon_deg * 10.0));
        for az_deg = 0:35
            fwd_az_rad = Nav_wrap_2pi(degtorad(az_deg * 10.0));
            for dist_m = -1:3
                dst_m = 10.0^dist_m;
                [ ep_lat_rad, ep_lon_rad, rev_az_rad, its ] = Nav_vincenty_direct(sp_lat_rad, sp_lon_rad, fwd_az_rad, dst_m, el, VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);
                [ adist_m, a_az_rad, a_raz_rad, its, is_ok ] = Nav_vincenty_inverse(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad, el, VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);
                
                if (abs(abs(a_az_rad - fwd_az_rad) - pi * 2) < 10E-3)
                    assert_approx_eq(abs(a_az_rad - fwd_az_rad), pi * 2, 10E-6);
                else
                    assert_approx_eq(fwd_az_rad, a_az_rad, 10E-6);
                end
                
                if abs(abs(rev_az_rad - a_raz_rad) - pi * 2) < 10E-3
                    assert_approx_eq(abs(rev_az_rad - a_raz_rad), pi * 2, 10E-6);           
                else
                    assert_approx_eq(rev_az_rad, a_raz_rad, 10E-6);
                end
                
                assert_approx_eq(dst_m, adist_m, 10E-6);
            end
        end
    end
end  
fprintf('Ok\n');


%%
fprintf('[ Testing ] Nav_test_toa_nlm_2d_solve without measurements noise...');
% random relevant values for actual point position
actual_x = 10.5;
actual_y = -11.0;
actual_z = 18.3;

x = -100.0;
y = 200.0;
z = 11.0;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z);
base_points(1,:) = [ x y z d ];

x = 115.0;
y = 267.0;
z = 17.3;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z);
base_points(2,:) = [ x y z d ];

x = 315.4;
y = -118.4;
z = 31.2;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z);
base_points(3,:) = [ x y z d ];

x = -470.1;
y = -216.7;
z = 12.5;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z);
base_points(4,:) = [ x y z d ];

[ x_best, y_best, rerr, it_cnt ] = Nav_toa_nlm_2d_solve(base_points, 0.0, 0.0, actual_z, NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 1.0);

assert_approx_eq(x_best, actual_x, 10E-4);
assert_approx_eq(y_best, actual_y, 10E-4);
assert_eq(rerr < 1.0, true, sprintf('Residual function value is greater than limit: %d > 1.0', rerr));
assert_eq(it_cnt < NLM_DEF_IT_LIMIT, true, sprintf('Number of iterations exceeded specified limit: %d >= %d', it_cnt, NLM_DEF_IT_LIMIT));
fprintf('Ok\n');

%%
fprintf('[ Testing ] Nav_test_toa_nlm_2d_solve with measurements noise...');
err_h_amp = 1.5; % distance measurement error in meters
% random relevant values for actual point position
actual_x = 10.5;
actual_y = -11.0;
actual_z = 18.3;

x = -100.0;
y = 200.0;
z = 11.0;
d_err = rand;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) + d_err * err_h_amp * 2.0 - err_h_amp;
base_points(1,:) = [ x y z d ];

x = 115.0;
y = 267.0;
z = 17.3;
d_err = rand;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) + d_err * err_h_amp * 2.0 - err_h_amp;
base_points(2,:) = [ x y z d ];

x = 315.4;
y = -118.4;
z = 31.2;
d_err = rand;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) + d_err * err_h_amp * 2.0 - err_h_amp;
base_points(3,:) = [ x y z d ];

x = -470.1;
y = -216.7;
z = 12.5;
d_err = rand;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) + d_err * err_h_amp * 2.0 - err_h_amp;
base_points(4,:) = [ x y z d ];

[ x_best, y_best, rerr, it_cnt ] = Nav_toa_nlm_2d_solve(base_points, 0.0, 0.0, actual_z, NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 1.0);

assert_approx_eq(x_best, actual_x, err_h_amp * 2.0);
assert_approx_eq(y_best, actual_y, err_h_amp * 2.0);
assert_eq(rerr < err_h_amp * 2.0, true, sprintf('Residual function value is greater than limit: %d > %d', rerr, err_h_amp * 2.0));
assert_eq(it_cnt < NLM_DEF_IT_LIMIT, true, sprintf('Number of iterations exceeded specified limit: %d >= %d', it_cnt, NLM_DEF_IT_LIMIT));
fprintf('Ok\n');

%%
fprintf('[ Testing ] Nav_test_tdoa_nlm_2d_solve with no measurements noise...');
% random relevant values for actual point position
actual_x = 10.5;
actual_y = -11.0;
actual_z = 18.3;               

base_points = zeros(4, 4);

x = -100.0;
y = 200.0;
z = 11.0;        
t = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) / propagation_velocity;
base_points(1,:) = [ x y z t ];

x = 115.0;
y = 267.0;
z = 17.3;        
t = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) / propagation_velocity;
base_points(2,:) = [ x y z t ];

x = 315.4;
y = -118.4;
z = 31.2;        
t = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) / propagation_velocity;
base_points(3,:) = [ x y z t ];

x = -470.1;
y = -216.7;
z = 12.5;        
t = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) / propagation_velocity;
base_points(4,:) = [ x y z t ];

base_lines = Nav_build_base_lines(base_points, propagation_velocity);

[ x_best, y_best, rerr, it_cnt ] = Nav_tdoa_nlm_2d_solve(base_lines, 0.0, 0.0, actual_z, NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 1.0);

assert_approx_eq(x_best, actual_x, 10E-4);
assert_approx_eq(y_best, actual_y, 10E-4);
assert_eq(rerr < 1.0, true, sprintf('Residual function value is greater than limit: %d > %d', rerr, 1.0));
assert_eq(it_cnt < NLM_DEF_IT_LIMIT, true, sprintf('Number of iterations exceeded specified limit: %d >= %d', it_cnt, NLM_DEF_IT_LIMIT));
fprintf('Ok\n');

%%
fprintf('[ Testing ] Nav_tdoa_nlm_3d_solve with no measurements noise...');
base_points = zeros(4, 4);
% random relevant values for actual point position
actual_x = 1;
actual_y = 2;
actual_z = 3;       

x = -60.0;
y = 20.0;
z = 10.0;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) / propagation_velocity;
base_points(1, :) = [ x y z d ];

x = 20.0;
y = -60.0;
z = 15.0;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) / propagation_velocity;
base_points(2, :) = [ x y z d ];

x = 20.0;
y = 60.0;
z = 20.0;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) / propagation_velocity;
base_points(3, :) = [ x y z d ];

x = -40.0;
y = -10.0;
z = 25.0;
d = Nav_dist_3d(x, y, z, actual_x, actual_y, actual_z) / propagation_velocity;
base_points(4, :) = [ x y z d ];

centroid = Nav_centroid(base_points);
base_lines = Nav_build_base_lines(base_points, propagation_velocity);

[ x_best, y_best, z_best, rerr, it_cnt ] = Nav_tdoa_nlm_3d_solve(base_lines, centroid(1), centroid(2), centroid(3), NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 10.0);

assert_approx_eq(x_best, actual_x, 1.0);
assert_approx_eq(y_best, actual_y, 1.0);
assert_approx_eq(z_best, actual_z, 5.0);

assert_eq(rerr < 1.0, true, sprintf('Residual function value is greater than limit: %d > %d', rerr, 1.0));
assert_eq(it_cnt < NLM_DEF_IT_LIMIT, true, sprintf('Number of iterations exceeded specified limit: %d >= %d', it_cnt, NLM_DEF_IT_LIMIT));
fprintf('Ok\n');


%% 
fprintf('[ Testing ] Nav_toa_locate_2d...');

base_number = 4;
start_base_z_m = 1.5;
base_z_step_m = 10.0;

actual_target_lat_deg = 44.505455; % singed degrees
actual_target_lon_deg = 48.547659; % signed degrees
actual_target_z_m = 25.0;          % meters

% generate base points via Vincenty equations
base_points = zeros(base_number, 4);

start_dst_projection_m = 500.0; % start distance from target, meters
dst_inc_step_m = 50.0;          % distance increment

% azimuth step
azimuth_step_rad = pi * 2 / base_number;

actual_target_lat_rad = degtorad(actual_target_lat_deg);
actual_target_lon_rad = degtorad(actual_target_lon_deg);

for base_idx = 1:base_number 

    dst_projection_m = start_dst_projection_m + dst_inc_step_m * (base_idx - 1);
    azimuth_rad = azimuth_step_rad * (base_idx - 1);

    [ep_lat_rad, ep_lon_rad, rev_az_rad, its] = Nav_vincenty_direct(actual_target_lat_rad, actual_target_lon_rad,...
        azimuth_rad, dst_projection_m,... 
        el,...
        VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);

    base_z_m = start_base_z_m + base_z_step_m * (base_idx - 1);
    dz_m = actual_target_z_m - base_z_m;    
    slant_range_m = sqrt(dst_projection_m * dst_projection_m + dz_m * dz_m);

    base_points(base_idx, :) = [ radtodeg(ep_lat_rad) radtodeg(ep_lon_rad) base_z_m slant_range_m ];
end

lat_prev_deg = nan;
lon_prev_deg = nan;

[ target_lat, target_lon, rerr, it_cnt ] = Nav_toa_locate_2d(base_points,...
    lat_prev_deg, lon_prev_deg, actual_target_z_m,...
    NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 10.0, el);

assert_approx_eq(actual_target_lat_deg, target_lat, 1E-6);
assert_approx_eq(actual_target_lon_deg, target_lon, 1E-6);
assert_eq(rerr < start_dst_projection_m * 0.01, true, sprintf('Residual function value is greater than limit: %d > %d', rerr, start_dst_projection_m * 0.01));
assert_eq(it_cnt < NLM_DEF_IT_LIMIT, true, sprintf('Number of iterations exceeded specified limit: %d >= %d', it_cnt, NLM_DEF_IT_LIMIT));

fprintf('Ok\n');

%% 
fprintf('[ Testing ] Nav_toa_locate_3d...');

base_number = 4;
start_base_z_m = 1.5;
base_z_step_m = 5.0;

actual_target_lat_deg = 44.505455; % singed degrees
actual_target_lon_deg = 48.547659; % signed degrees
actual_target_z_m = 25.0;          % meters

% generate base points via Vincenty equations
base_points = zeros(base_number, 4);

start_dst_projection_m = 500.0; % start distance from target, meters
dst_inc_step_m = 50.0;          % distance increment

% azimuth step
azimuth_step_rad = 2 * pi / base_number;

actual_target_lat_rad = degtorad(actual_target_lat_deg);
actual_target_lon_rad = degtorad(actual_target_lon_deg);        

for base_idx = 1:base_number
    dst_projection_m = start_dst_projection_m + dst_inc_step_m * (base_idx - 1);
    azimuth_rad = azimuth_step_rad * (base_idx - 1);

    [ep_lat_rad, ep_lon_rad, rev_az_rad, its] = Nav_vincenty_direct(actual_target_lat_rad, actual_target_lon_rad,...
        azimuth_rad, dst_projection_m,...
        el,...
        VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);

    base_z_m = start_base_z_m + base_z_step_m * (base_idx - 1);
    dz_m = actual_target_z_m - base_z_m;    
    slant_range_m = sqrt(dst_projection_m * dst_projection_m + dz_m * dz_m);

    base_points(base_idx, :) = [ radtodeg(ep_lat_rad) radtodeg(ep_lon_rad) base_z_m slant_range_m ];
end

lat_prev_deg = nan;
lon_prev_deg = nan;
z_prev_m = nan;
[ target_lat_deg, target_lon_deg, target_z_m, rerr, it_cnt ] = Nav_toa_locate_3d(base_points,...
    lat_prev_deg, lon_prev_deg, z_prev_m,...
    NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 10.0, el);

assert_approx_eq(actual_target_lat_deg, target_lat_deg, 1E-6);
assert_approx_eq(actual_target_lon_deg, target_lon_deg, 1E-6);
assert_eq(rerr < start_dst_projection_m * 0.01, true, sprintf('Residual function value is greater than limit: %d > %d', rerr, start_dst_projection_m * 0.01));
assert_eq(it_cnt < NLM_DEF_IT_LIMIT, true, sprintf('Number of iterations exceeded specified limit: %d >= %d', it_cnt, NLM_DEF_IT_LIMIT));

fprintf('Ok\n');

%%
fprintf('[ Testing ] Nav_tdoa_locate_2d...');

base_number = 4;
start_base_z_m = 1.5;
base_z_step_m = 10.0;

actual_target_lat_deg = 44.505455; % singed degrees
actual_target_lon_deg = 48.547659; % signed degrees
actual_target_z_m = 25.0;          % meters

% generate base points via Vincenty equations
base_points = zeros(base_number, 4);

start_dst_projection_m = 500.0; % start distance from target, meters
dst_inc_step_m = 50.0;          % distance increment

% azimuth step
azimuth_step_rad = pi * 2 / base_number;

actual_target_lat_rad = degtorad(actual_target_lat_deg);
actual_target_lon_rad = degtorad(actual_target_lon_deg);

% signal propagation speed

for base_idx = 1:base_number

    dst_projection_m = start_dst_projection_m + dst_inc_step_m * (base_idx - 1);
    azimuth_rad = azimuth_step_rad * (base_idx - 1);

    [ ep_lat_rad, ep_lon_rad, rev_az_rad, its ] = Nav_vincenty_direct(actual_target_lat_rad, actual_target_lon_rad,... 
        azimuth_rad, dst_projection_m,...
        el,...
        VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);

    base_z_m = start_base_z_m + base_z_step_m * (base_idx - 1);
    dz_m = actual_target_z_m - base_z_m;    
    slant_range_m = sqrt(dst_projection_m * dst_projection_m + dz_m * dz_m);

    base_points(base_idx, :) = [ radtodeg(ep_lat_rad) radtodeg(ep_lon_rad) base_z_m slant_range_m / propagation_velocity ];
end


lat_prev_deg = nan;
lon_prev_deg = nan;
[ target_lat_deg, target_lon_deg, rerr, it_cnt ] = Nav_tdoa_locate_2d(base_points,...
    lat_prev_deg, lon_prev_deg, actual_target_z_m,...
    NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 10.0, el, propagation_velocity);

assert_approx_eq(actual_target_lat_deg, target_lat_deg, 1E-6);
assert_approx_eq(actual_target_lon_deg,target_lon_deg, 1E-6);
assert_eq(rerr < start_dst_projection_m * 0.01, true, sprintf('Residual function value is greater than limit: %d > %d', rerr, start_dst_projection_m * 0.01));
assert_eq(it_cnt < NLM_DEF_IT_LIMIT, true, sprintf('Number of iterations exceeded specified limit: %d >= %d', it_cnt, NLM_DEF_IT_LIMIT));

fprintf('Ok\n');


%%
fprintf('[ Testing ] Nav_tdoa_locate_3d...');

base_number = 4;
start_base_z_m = 1.5;
base_z_step_m = 5.0;

actual_target_lat_deg = 48.513724; % singed degrees
actual_target_lon_deg = 44.553248; % signed degrees
actual_target_z_m = 25.0;          % meters

% generate base points via Vincenty equations
base_points = zeros(base_number, 4);

start_dst_projection_m = 500.0;          % start distance from target, meters
dst_inc_step_m = 50.0;                   % distance increment

% azimuth step
azimuth_step_rad = 2 * pi / base_number;

actual_target_lat_rad = degtorad(actual_target_lat_deg);
actual_target_lon_rad = degtorad(actual_target_lon_deg);

for base_idx = 1:base_number

    dst_projection_m = start_dst_projection_m + dst_inc_step_m * (base_idx - 1);
    azimuth_rad = azimuth_step_rad * (base_idx - 1);

    [ ep_lat_rad, ep_lon_rad, rev_az_rad, its ] = Nav_vincenty_direct(actual_target_lat_rad, actual_target_lon_rad,...
        azimuth_rad, dst_projection_m,...
        el,...
        VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);

    base_z_m = start_base_z_m + base_z_step_m * (base_idx - 1);
    dz_m = actual_target_z_m - base_z_m;    
    slant_range_m = sqrt(dst_projection_m * dst_projection_m + dz_m * dz_m);

    base_points(base_idx, :) = [ radtodeg(ep_lat_rad) radtodeg(ep_lon_rad) base_z_m (slant_range_m / propagation_velocity) ];
end


lat_prev_deg = nan;
lon_prev_deg = nan;
prev_z_m = nan;

[ target_lat_deg, target_lon_deg, target_z_m, rerr, it_cnt ] = Nav_tdoa_locate_3d(base_points,...
    lat_prev_deg, lon_prev_deg, prev_z_m,...
    NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 10.0, el, propagation_velocity);

[ dst_m, fwd_az_rad, rev_az_rad, its, is_ok ] = Nav_vincenty_inverse(actual_target_lat_rad, actual_target_lon_rad,... 
    target_lat_deg * pi / 180, target_lon_deg * pi / 180,...
    el, VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);

assert_eq(dst_m < start_dst_projection_m * 0.01, true, sprintf('Estimated location is further than limit (%d): %d', start_dst_projection_m * 0.01, dst_m));
assert_approx_eq(target_z_m, actual_target_z_m, start_dst_projection_m * 0.05);
assert_eq(rerr < start_dst_projection_m * 0.01, true, sprintf('Residual function value is greater than limit: %d > %d', rerr, start_dst_projection_m * 0.01));
assert_eq(it_cnt < NLM_DEF_IT_LIMIT, true, sprintf('Number of iterations exceeded specified limit: %d >= %d', it_cnt, NLM_DEF_IT_LIMIT));

fprintf('Ok\n');




%% 
fprintf('[ Testing ] Nav_aoa_ng_2d_solve...');

aoa_target_range_m = 1000;
aoa_n_base_points = 4;
aoa_base_size_m = 3;
aoa_base_points = zeros(aoa_n_base_points, 4);
aoa_deg_error = 0.1;

for n = 1:aoa_n_base_points   
   az_ = (n - 1) * 2 * pi / (aoa_n_base_points);
   base_x = aoa_base_size_m * cos(az_);
   base_y = aoa_base_size_m * sin(az_);
   base_z = 0;
   aoa_base_points(n, :) = [ base_x base_y base_z 0 ];   
end

for actual_angle = 0:359
   actual_angle_rad = actual_angle * pi / 180;
   target_x = aoa_target_range_m * cos(actual_angle_rad);
   target_y = aoa_target_range_m * sin(actual_angle_rad);
   target_z = 0;
          
   for n = 1:aoa_n_base_points
       aoa_base_points(n, 4) = Nav_dist_3d(aoa_base_points(n, 1), aoa_base_points(n, 2), aoa_base_points(n, 3),...
           target_x, target_y, target_z) / propagation_velocity;              
   end
      
   [a_rad, it_cnt] = Nav_tdoa_aoa_ng_2d_solve(aoa_base_points, AOA_DEF_IT_LIMIT,...
     AOA_DEF_PREC_THRLD, propagation_velocity);
   
   aoa_angular_error = actual_angle_rad - a_rad;
   if (aoa_angular_error > pi)
      aoa_angular_error = 2 * pi - aoa_angular_error; 
   end
   
   aoa_angular_error = aoa_angular_error * 180 / pi;
 
   assert_eq(it_cnt < AOA_DEF_IT_LIMIT, true, sprintf('Number of iterations exceeded specified limit: %d >= %d', it_cnt, AOA_DEF_IT_LIMIT));
   assert_eq(aoa_angular_error < aoa_deg_error, true, sprintf('AOA estimation error exeeded specified limit: %d >= %d', aoa_angular_error, aoa_deg_error));
end

fprintf('Ok\n');


%%
timeElapsed = toc;
fprintf('\nAll tests passed successfully!\n');
fprintf('Elapsed time: %.3f seconds\n', timeElapsed);



