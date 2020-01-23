function [ new_lat_rad, new_lon_rad ] = Nav_geopoint_offset_by_deltas_wgs84(lat_rad, lon_rad, lat_offset_m, lon_offset_m)

    m_per_deg_lat = 111132.92 - 559.82 * cos(2.0 * lat_rad) + 1.175 * cos(4.0 * lat_rad);
    m_per_deg_lon = 111412.84 * cos(lat_rad) - 93.5 * cos(3.0 * lat_rad);
    
    new_lat_rad = lat_rad - (pi / 180) * lat_offset_m / m_per_deg_lat;
    new_lon_rad = lon_rad - (pi / 180) * lon_offset_m / m_per_deg_lon;
end