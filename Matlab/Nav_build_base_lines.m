function [ base_lines ] = Nav_build_base_lines(base_points, velocity)

    [n_base_points, nc] = size(base_points);
    base_lines = zeros(1, 7);

    k_ = 1;
    for i_ = 1:n_base_points - 1
        for j_ = i_ + 1:n_base_points
            
            base_lines(k_,:) = [ base_points(i_, 1) base_points(i_, 2) base_points(i_, 3)...
                                base_points(j_, 1) base_points(j_, 2) base_points(j_, 3)...
                                (base_points(i_, 4) - base_points(j_, 4)) * velocity ];
                            
            k_ = k_ + 1;
        end
    end  

end