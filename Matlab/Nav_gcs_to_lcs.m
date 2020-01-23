% Converts geographic coordinates to local cartesian coordinates
function [ metric_points ] = Nav_gcs_to_lcs(points, el)

[ n_points, n_coords ] = size(points);
metric_points = zeros(n_points, n_coords);

centroid_deg = Nav_centroid(points);
centroid_rad = centroid_deg * (pi / 180);

for n = 1:n_points
    [ d_lat_m, d_lon_m ] = Nav_get_deltas_by_geopoints(centroid_rad(1), centroid_rad(2),...
        points(n, 1) * pi / 180, points(n, 2) * pi / 180, el);
    metric_points(n, :) = [ d_lat_m d_lon_m points(n, 3:end) ];
end

end