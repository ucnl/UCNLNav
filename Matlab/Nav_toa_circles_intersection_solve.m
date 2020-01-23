function [ a_best ] = Nav_toa_circles_intersection_solve(base_points, anchor_x, anchor_y, radius, z,...
                            arc_mid_rad, arc_angle_rad, steps)

a = arc_mid_rad - arc_angle_rad / 2.0;
a_end = arc_mid_rad + arc_angle_rad / 2.0;
a_best = a;
step_rad = arc_angle_rad / steps;
eps_best = realmax;

while a < a_end
    x = anchor_x + radius * cos(a);
    y = anchor_y + radius * sin(a);
    eps = Nav_eps_toa3d(base_points, x, y, z);

    if eps < eps_best
        eps_best = eps;
        a_best = a;
    end

    a = a + step_rad;
end

end