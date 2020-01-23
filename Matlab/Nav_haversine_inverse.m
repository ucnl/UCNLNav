function [ result ] = Nav_haversine_inverse(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad, e_radius_m)
    
   a = sin((ep_lat_rad - sp_lat_rad) / 2.0)^2 + cos(sp_lat_rad) * cos(ep_lat_rad) *...
       sin((ep_lon_rad - sp_lon_rad) / 2.0)^2;
            
   result = e_radius_m * 2.0 * atan2( sqrt(a), sqrt(1.0 - a));
end