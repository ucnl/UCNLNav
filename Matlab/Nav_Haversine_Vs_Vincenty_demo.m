close all;
clear all;
clc;

% Nav_Haversine_Vs_Vincenty_demo.m

VNC_DEF_EPSILON = 1E-12;
VNC_DEF_IT_LIMIT = 2000;

sp_lat_rad = degtorad(48.527683);
sp_lon_rad = degtorad(44.558815);
fwd_az_rad = 1.5 * pi + (rand * pi / 4 - pi / 8);

el = Nav_build_standard_ellipsoid('WGS84');


n_samples = 10000;
v_lats_rad = zeros(n_samples, 1);
v_lons_rad = zeros(n_samples, 1);
h_lats_rad = zeros(n_samples, 1);
h_lons_rad = zeros(n_samples, 1);
v_dist = zeros(n_samples, 1);
ip_v_dist = zeros(n_samples, 1);
ip_h_dist = zeros(n_samples, 1);

%% For long distance (up to 10000 km)

step_m = 1000; % meters
distances = (1:n_samples) .* step_m;

for idx = 1:n_samples
           
        [ h_lats_rad(idx), h_lons_rad(idx) ] = Nav_haversine_direct(sp_lat_rad,...
            sp_lon_rad,...
            distances(idx),...
            fwd_az_rad,...
            el.mjsa_m);
        
        [ v_lats_rad(idx), v_lons_rad(idx), v_rev_az_rad, v_its ] = Nav_vincenty_direct(sp_lat_rad,...
            sp_lon_rad,...
            fwd_az_rad,...
            distances(idx),...
            el,...
            VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);
        
        [ v_dist(idx) a_az_rad, a_raz_rad, its, is_ok ] = Nav_vincenty_inverse(h_lats_rad(idx),...
            h_lons_rad(idx),...
            v_lats_rad(idx),...
            v_lons_rad(idx),...
            el,...
            VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);
        
        [ ip_v_dist(idx) a_az_rad, a_raz_rad, its, is_ok ] = Nav_vincenty_inverse(sp_lat_rad,...
            sp_lon_rad,...
            v_lats_rad(idx),...
            v_lons_rad(idx),...
            el,...
            VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);
        
        ip_h_dist(idx) = Nav_haversine_inverse(sp_lat_rad,...
            sp_lon_rad,...
            v_lats_rad(idx),...
            v_lons_rad(idx),...
            el.mjsa_m);
end


figure
hold on
[vx, vy, vz] = sph2cart(v_lons_rad, v_lats_rad, el.mjsa_m);
plot3(vx, vy, vz, '-r');
[hx, hy, hz] = sph2cart(h_lons_rad, h_lats_rad, el.mjsa_m);
plot3(hx, hy, hz, '-b');

latspacing = 20; 
lonspacing = 20;


[lon1, lat1] = meshgrid(-180:lonspacing: 180,linspace(-90,90,300)); 
[x1, y1, z1] = sph2cart(lon1 * pi / 180, lat1 * pi / 180, el.mjsa_m); 
plot3(x1,y1,z1,'-', 'color', [ 0.5 0.5 0.5 ]);

[lat2, lon2] = meshgrid(-90:latspacing:90, linspace(-180,180,300)); 
[x2, y2, z2] = sph2cart(lon2 * pi / 180, lat2 * pi / 180, el.mjsa_m); 
plot3(x2,y2,z2,'-', 'color', [ 0.5 0.5 0.5 ]);
axis equal tight off
title('Direct geodetic problem: Vincenty vs. Haversine');
legend('Vincenty', 'Haversine');
view(90, 22);


d_lat_deg = radtodeg(v_lats_rad - h_lats_rad);
d_lon_deg = radtodeg(v_lons_rad - h_lons_rad);



figure
semilogx(distances, d_lat_deg, 'r');
title('Direct geodetic problem: Vincenty vs. Haversine (Latitude difference)');
xlabel('Distance, m');
ylabel('Difference, °');

figure
semilogx(distances, d_lon_deg, 'r');
title('Direct geodetic problem: Vincenty vs. Haversine (Longitude difference)');
xlabel('Distance, m');
ylabel('Difference, °');

figure
semilogx(distances, v_dist, 'r');
title('Direct geodetic problem: Vincenty vs. Haversine (Endpoint difference by Vincenty)');
xlabel('Distance, m');
ylabel('Difference, m');

figure
semilogx(distances, ip_v_dist - ip_h_dist, 'r');
title('Inverse geodetic problem: Vincenty vs. Haversine (Distance difference)');
xlabel('Distance, m');
ylabel('Difference, m');



%% For small distances (up to 10 km)

step_m = 1; % meters
distances = (1:n_samples) .* step_m;

for idx = 1:n_samples
           
        [ h_lats_rad(idx), h_lons_rad(idx) ] = Nav_haversine_direct(sp_lat_rad,...
            sp_lon_rad,...
            distances(idx),...
            fwd_az_rad,...
            el.mjsa_m);
        
        [ v_lats_rad(idx), v_lons_rad(idx), v_rev_az_rad, v_its ] = Nav_vincenty_direct(sp_lat_rad,...
            sp_lon_rad,...
            fwd_az_rad,...
            distances(idx),...
            el,...
            VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);
        
        [ v_dist(idx) a_az_rad, a_raz_rad, its, is_ok ] = Nav_vincenty_inverse(h_lats_rad(idx),...
            h_lons_rad(idx),...
            v_lats_rad(idx),...
            v_lons_rad(idx),...
            el,...
            VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);
        
        [ ip_v_dist(idx) a_az_rad, a_raz_rad, its, is_ok ] = Nav_vincenty_inverse(sp_lat_rad,...
            sp_lon_rad,...
            v_lats_rad(idx),...
            v_lons_rad(idx),...
            el,...
            VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);
        
        ip_h_dist(idx) = Nav_haversine_inverse(sp_lat_rad,...
            sp_lon_rad,...
            v_lats_rad(idx),...
            v_lons_rad(idx),...
            el.mjsa_m);
end


d_lat_deg = radtodeg(v_lats_rad - h_lats_rad);
d_lon_deg = radtodeg(v_lons_rad - h_lons_rad);

figure
semilogx(distances, d_lat_deg, 'r');
title('Direct geodetic problem: Vincenty vs. Haversine (Latitude difference)');
xlabel('Distance, m');
ylabel('Difference, °');

figure
semilogx(distances, d_lon_deg, 'r');
title('Direct geodetic problem: Vincenty vs. Haversine (Longitude difference)');
xlabel('Distance, m');
ylabel('Difference, °');

figure
semilogx(distances, v_dist, 'r');
title('Direct geodetic problem: Vincenty vs. Haversine (Endpoint difference by Vincenty)');
xlabel('Distance, m');
ylabel('Difference, m');

figure
semilogx(distances, ip_v_dist - ip_h_dist, 'r');
title('Inverse geodetic problem: Vincenty vs. Haversine (Distance difference)');
xlabel('Distance, m');
ylabel('Difference, m');



