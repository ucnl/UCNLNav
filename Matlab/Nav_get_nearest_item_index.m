
function [ nrst_idx ] = Nav_get_nearest_item_index(base_points)

[ n_points, n_coords ] = size(base_points);
nrst_idx = 1;
min_dst = base_points(1, 4);

for idx = 1:n_points
    if base_points(idx, 4) < min_dst 
        min_dst = base_points(idx, 4);
        nrst_idx = idx;
    end
end

end