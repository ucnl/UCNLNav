function [ result ] = Nav_dist_3d(x1, y1, z1, x2, y2, z2)
   
result = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);

end