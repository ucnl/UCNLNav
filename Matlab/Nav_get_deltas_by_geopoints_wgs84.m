function [ delta_lat_m, delta_lon_m ] = Nav_get_deltas_by_geopoints_wgs84(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad)   
    
    m_lat_rad = (sp_lat_rad + ep_lat_rad) / 2.0;
    m_per_deg_lat = 111132.92 - 559.82 * cos(2.0 * m_lat_rad) + 1.175 * cos(4.0 * m_lat_rad);
    m_per_deg_lon = 111412.84 * cos(m_lat_rad) - 93.5 * cos(3.0 * m_lat_rad);
    
    delta_lat_m = (sp_lat_rad - ep_lat_rad) * m_per_deg_lat * (180 / pi);
    delta_lon_m = (sp_lon_rad - ep_lon_rad) * m_per_deg_lon * (180 / pi);
end