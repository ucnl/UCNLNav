function [ result ] = Nav_haversine_final_bearing(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad)
   result = Nav_wrap_2pi(pi + Nav_haversine_initial_bearing(ep_lat_rad, ep_lon_rad, sp_lat_rad, sp_lon_rad));
end
