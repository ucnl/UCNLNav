function [ new_lat_rad, new_lon_rad ] = Nav_geopoint_offset_by_deltas(lat_rad, lon_rad, lat_offset_m, lon_offset_m, el)
    m_per_deg_lat = Nav_lat_1deg_length(lat_rad, el);
    m_per_deg_lon = Nav_lon_1deg_length(lat_rad, el);
    
    new_lat_rad = lat_rad - (pi / 180) * lat_offset_m / m_per_deg_lat;
    new_lon_rad = lon_rad - (pi / 180) * lon_offset_m / m_per_deg_lon;
end 