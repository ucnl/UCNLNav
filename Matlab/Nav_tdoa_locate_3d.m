% Solves TDOA navigation problem: searches for a target location by base points - points with known locations
% and measured times of arrival. Nelder-Mead (simplex) method is used.
function [ target_lat_deg, target_lon_deg, target_z_m, rerr, it_cnt ] = Nav_tdoa_locate_3d(bases,...
            lat_prev_deg, lon_prev_deg, z_prev_m,...
            max_iterations, precision_threshold, simplex_size,...
            el, velocity)

bases_centroid_deg = Nav_centroid(bases);
bases_centroid_rad = bases_centroid_deg * pi / 180;
base_points_m = Nav_gcs_to_lcs(bases, el);
base_lines_m = Nav_build_base_lines(base_points_m, velocity);

x_prev_m = 0.0;
y_prev_m = 0.0;
z_prev_m_ = z_prev_m;

if isnan(z_prev_m)
    z_prev_m_ = bases_centroid_deg(3);
end

if ~isnan(lat_prev_deg) && ~isnan(lon_prev_deg)
    [ x_prev_m, y_prev_m ] = Nav_get_deltas_by_geopoints(bases_centroid_rad(1), bases_centroid_rad(2),... 
        lat_prev_deg * pi / 180, lon_prev_deg * pi / 180, el);
 end

[ x_best, y_best, target_z_m, rerr, it_cnt ] = Nav_tdoa_nlm_3d_solve(base_lines_m, x_prev_m, y_prev_m, z_prev_m_,...
    max_iterations, precision_threshold, simplex_size);

[ target_lat_rad, target_lon_rad ] = Nav_geopoint_offset_by_deltas(bases_centroid_rad(1), bases_centroid_rad(2),...
    y_best, x_best, el);

target_lat_deg = target_lat_rad * 180 / pi;
target_lon_deg = target_lon_rad * 180 / pi;

end