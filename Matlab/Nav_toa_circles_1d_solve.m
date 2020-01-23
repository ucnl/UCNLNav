function [ x_best, y_best, rerr ] = Nav_toa_circles_1d_solve(base_points, z,...
    end_arc_angle_rad, steps, arc_angle_decrease_factor)

nrst_idx = Nav_get_nearest_item_index(base_points);
d_z = abs(base_points(nrst_idx, 3) - z);

if base_points(nrst_idx, 4) < d_z 
    radius = 0;
else
    radius = sqrt(base_points(nrst_idx, 4)^2 - d_z^2);
end

anchor_x = base_points(nrst_idx, 1);
anchor_y = base_points(nrst_idx, 2);
arc_angle = pi * 2;
alpha = 0.0;    

while arc_angle > end_arc_angle_rad 
    alpha = Nav_toa_circles_intersection_solve(base_points, anchor_x, anchor_y, radius, z, alpha, arc_angle, steps);
    arc_angle = arc_angle * arc_angle_decrease_factor;
end

x_best = anchor_x + radius * cos(alpha);
y_best = anchor_y + radius * sin(alpha);

rerr = sqrt(Nav_eps_toa3d(base_points, x_best, y_best, z));

end