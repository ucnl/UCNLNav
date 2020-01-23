function [ delta_lat_m, delta_lon_m ] = Nav_get_deltas_by_geopoints(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad, el)
    
    m_lat_rad = (sp_lat_rad + ep_lat_rad) / 2.0;
    m_per_deg_lat = Nav_lat_1deg_length(m_lat_rad, el);
    m_per_deg_lon = Nav_lon_1deg_length(m_lat_rad, el);
    
    delta_lat_m = (sp_lat_rad - ep_lat_rad) * m_per_deg_lat * (180 / pi);
    delta_lon_m = (sp_lon_rad - ep_lon_rad) * m_per_deg_lon * (180 / pi);
end