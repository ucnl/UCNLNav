function [ dst_m, fwd_az_rad, rev_az_rad, its, is_ok ] = Nav_vincenty_inverse(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad, el, epsilon, it_limit) %-> (f64, f64, f64, i32, bool)
    
    l_ = ep_lon_rad - sp_lon_rad;
    tan_u_1 = (1.0 - el.fltn) * tan(sp_lat_rad);
    cos_u_1 = 1.0 / sqrt(1.0 + tan_u_1^2);
    sin_u_1 = tan_u_1 * cos_u_1;
    tan_u_2 = (1.0 - el.fltn) * tan(ep_lat_rad);
    cos_u_2 = 1.0 / sqrt(1.0 + tan_u_2^2);
    sin_u_2 = tan_u_2 * cos_u_2;
            
    sin_sigma = 0.0;
    cos_sigma = 0.0;
    cos_sq_alpha = 0.0;
    cos_2_sigma_m = 0.0;
    sigma = 0.0;

    lambda = l_;
    lambda_ = 0.0;
    its = 0;

    it_check = 0.0;    
    antimeridian = abs(l_) > pi;

    do = 1;
    
    while do

        sin_lambda = sin(lambda);
        cos_lambda = cos(lambda);

        sin_sq_sigma = (cos_u_2 * sin_lambda) * (cos_u_2 * sin_lambda) +...
                       (cos_u_1 * sin_u_2 - sin_u_1 * cos_u_2 * cos_lambda) *...
                       (cos_u_1 * sin_u_2 - sin_u_1 * cos_u_2 * cos_lambda);

        if abs(sin_sq_sigma) > eps
        
            sin_sigma = sqrt(sin_sq_sigma);
            cos_sigma = sin_u_1 * sin_u_2 + cos_u_1 * cos_u_2 * cos_lambda;
            sigma = atan2(sin_sigma, cos_sigma);
            sin_alpha = cos_u_1 * cos_u_2 * sin_lambda / sin_sigma;

            cos_sq_alpha = 1.0 - sin_alpha^2;
            
            if (cos_sq_alpha ~= 0.0)
                cos_2_sigma_m = cos_sigma - 2.0 * sin_u_1 * sin_u_2 / cos_sq_alpha;            
            else
                cos_sq_alpha = 0.0;
            end

            c_ = el.fltn / 16.0 * cos_sq_alpha * (4.0 + el.fltn * (4.0 - 3.0 * cos_sq_alpha));
            lambda_ = lambda;
            lambda = l_ + (1.0 - c_) * el.fltn * sin_alpha *...
                     (sigma + c_ * sin_sigma * (cos_2_sigma_m + c_ * cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m * cos_2_sigma_m)));
        
            if (antimeridian)
                it_check = abs(lambda) - pi;
            else
                it_check = abs(lambda);
            end     
        end
            
        its = its + 1;
        do = ((abs(lambda - lambda_) > epsilon) && (its < it_limit) && (it_check < pi));
    end

    u_sq = cos_sq_alpha * (el.mjsa_m^2 - el.mnsa_m^2) / el.mnsa_m^2;
    a_ = 1.0 + u_sq / 16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
    b_ = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));
    delta_sigma = b_ * sin_sigma * (cos_2_sigma_m + b_/4.0 * (cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m * cos_2_sigma_m) -...
        b_/6.0 * cos_2_sigma_m * (-3.0 + 4.0 * sin_sigma * sin_sigma) * (-3.0 + 4.0 * cos_2_sigma_m * cos_2_sigma_m)));

    dst_m = el.mnsa_m * a_ * (sigma - delta_sigma);

    fwd_az_rad = Nav_wrap_2pi(atan2(cos_u_2 * sin_lambda, cos_u_1 * sin_u_2 - sin_u_1 * cos_u_2 * cos_lambda));
    rev_az_rad = Nav_wrap_2pi(atan2(cos_u_1 * sin_lambda, - sin_u_1 * cos_u_2 + cos_u_1 * sin_u_2 * cos_lambda));
    
    is_ok = (its < it_limit) && (it_check < pi);

end

