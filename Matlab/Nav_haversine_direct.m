function [ ep_lat_rad, ep_lon_rad ] = Nav_haversine_direct(sp_lat_rad, sp_lon_rad, dst_m, fwd_az_rad, e_radius_m)

    delta = dst_m / e_radius_m;
    
    ep_lat_rad = Nav_wrap_2pi(asin(sin(sp_lat_rad) * cos(delta) + cos(sp_lat_rad) * sin(delta) * cos(fwd_az_rad)));
    
    ep_lon_rad = Nav_wrap_2pi(3 * pi +... 
        (sp_lon_rad + atan2(sin(fwd_az_rad) * sin(delta) * cos(sp_lat_rad), (cos(delta) - sin(sp_lat_rad) * sin(ep_lat_rad))))) - pi;
end
