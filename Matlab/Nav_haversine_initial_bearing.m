        
function [ result ] = Nav_haversine_initial_bearing(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad)
    
    y = sin(ep_lon_rad - sp_lon_rad) * cos(ep_lat_rad);
    x = cos(sp_lat_rad) * sin(ep_lat_rad) - sin(sp_lat_rad) * cos(ep_lat_rad) * cos(ep_lon_rad - sp_lon_rad);
    
    result = Nav_wrap_2pi(atan2(y, x));
end