% base_points(n, c)
% n - a base point index
% c = 1 -> x
% c = 2 -> y
% c = 3 -> z
% c = 4 -> estimated distance

function [ result ] = Nav_eps_toa3d(base_points, x, y, z)
    
result = 0.0;
    
for n = 1:length(base_points)
   
    result = result + (sqrt((base_points(n, 1) - x)^2 +...
                            (base_points(n, 2) - y)^2 +...
                            (base_points(n, 3) - z)^2) - base_points(n, 4))^2;
end