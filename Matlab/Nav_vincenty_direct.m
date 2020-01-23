function [ ep_lat_rad, ep_lon_rad, rev_az_rad, its ] = Nav_vincenty_direct(sp_lat_rad, sp_lon_rad, fwd_az_rad, dst_m, el, epsilon, it_limit)

    sin_alpha_1 = sin(fwd_az_rad);
    cos_alpha_1 = cos(fwd_az_rad);
    tan_u_1 = (1.0 - el.fltn) * tan(sp_lat_rad);
    cos_u_1 = 1.0 / sqrt(1.0 + tan_u_1^2);
    sin_u_1 = tan_u_1 * cos_u_1;

    sigma_1 = atan2(tan_u_1, cos_alpha_1);
    sin_alpha = cos_u_1 * sin_alpha_1;
    cos_sq_alpha = 1.0 - sin_alpha^2;
    u_sq = cos_sq_alpha * (el.mjsa_m^2 - el.mnsa_m^2) / el.mnsa_m^2;
    a_ = 1.0 + u_sq / 16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
    b_ = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));
    
    sigma = dst_m / (el.mnsa_m * a_);
    its = 0;
            
    do = 1;
    
    while do
        
        cos_2_sigma_m = cos(2.0 * sigma_1 + sigma);
        sin_sigma = sin(sigma);
        cos_sigma = cos(sigma);

        delta_sigma = b_ * sin_sigma * (cos_2_sigma_m + b_ / 4.0 * (cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m * cos_2_sigma_m) -...
                      b_ / 6.0 * cos_2_sigma_m * (-3.0 + 4.0 * sin_sigma * sin_sigma) * (-3.0 + 4.0 * cos_2_sigma_m * cos_2_sigma_m)));

        sigma_ = sigma;
        sigma = dst_m / (el.mnsa_m * a_) + delta_sigma;

        its = its + 1;
        do = ((abs(sigma - sigma_) > epsilon) && (its < it_limit));
    end
    
    x = sin_u_1 * sin_sigma - cos_u_1 * cos_sigma * cos_alpha_1;
    
    ep_lat_rad = Nav_wrap_2pi(atan2(sin_u_1 * cos_sigma + cos_u_1 * sin_sigma * cos_alpha_1, (1.0 - el.fltn) * sqrt(sin_alpha^2 + x^2)));
    
    lambda = atan2(sin_sigma * sin_alpha_1, cos_u_1 * cos_sigma - sin_u_1 * sin_sigma * cos_alpha_1);
    c_ = el.fltn / 16.0 * cos_sq_alpha * (4.0 + el.fltn * (4.0 - 3.0 * cos_sq_alpha));    
    l_ = lambda - (1.0 - c_) * el.fltn * sin_alpha * (sigma + c_ * sin_sigma *... 
             (cos_2_sigma_m + c_ * cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m^2)));
    
    ep_lon_rad = Nav_wrap_2pi(sp_lon_rad + l_);
    rev_az_rad = Nav_wrap_2pi(atan2(sin_alpha, -x));

end
