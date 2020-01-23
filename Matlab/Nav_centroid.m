function [ centroid ] = Nav_centroid(points)
    
[ n_points, n_coords ] = size(points);
centroid = zeros(n_coords, 1);

for i_=1:n_points
    for j_ = 1:n_coords
        centroid(j_) = centroid(j_) + points(i_, j_);
    end
end

centroid = centroid ./ n_points;

end