% Solves the AOA navigation problem: searches for a horizontal angle of
% arrival by base points - points with known locations and measured times
% of arrival. Newton-Gauss method is used.
 function [a_rad, it_cnt] = Nav_tdoa_aoa_ng_2d_solve(base_points, max_iterations, ...
     precision_threshold, velocity)

% Source function:
% eps(a,e) = Sum[dXij*cose*cosa + dYij*cose*sina + dZij*sine - dTij*v]^2
%
% functional:
% fij(a,e) = dXij*cose*cosa + dYij*cose*sina + dZij*sine - dTij*v
%
% Linearized by Taylor series decomposition
% fij(a,e) ~ fij(a0,e0) + dF(a,e)/da * (a - a0) + dF(a,e)/de * (e - e0) = 0
% 
% dF(a,e)/da = dXij*cose*(-sina) + dYij*cose*cosa
% dF(a,e)/de = dXij(-sine)*cosa + dYij*(-sine)sina + dZij*cose
%
%
% Y = | a |
%     | e |
%
% X = | dFij(a,e)/da     dFij(a,e)/de |
%     | ...              ...          |
%     | ...              ...          |
%
% b = | fij(a0,e0) |
%     | ...        |
%     | ...        |
%
% Y = (X^T * X)^-1 * X^T * b
%
%
%
% For horizontal angle only:
%
% Y = | a |
%
% X = | dFij(a)/da |
%     | ...        |
%     | ...        |
%
% b = | fij(a0)    |
%     | ...        |
%     | ...        |
%
% Y = X^T/(X^T * X) * b
% 

% estimate the earliest sensor
[n_sensors, nc] = size(base_points);

earliest_sensor_idx = 1;
earliest_toa = base_points(earliest_sensor_idx, 4);
for n = 2:n_sensors
    if (base_points(n, 4) < earliest_toa)
       earliest_toa = base_points(n, 4);
       earliest_sensor_idx = n;
    end
end

% the number of time-differencies
dt_number = factorial(n_sensors) / ((n_sensors - 2) * factorial(2));

% estimate the first approximation
%
%             a = 90
%             |
%     i=1     |       i+1
%             |
%             |
%             +------------- a = 0
%             |
%             |
%     i+...   |       i+2
%             |
%
a0 = atan2(base_points(earliest_sensor_idx, 2), base_points(earliest_sensor_idx, 1));

a_rad = a0;
it_cnt = 0;

dfda = zeros(dt_number, 1);
dt = zeros(dt_number, 1);
dr = zeros(dt_number, 1);
dx = zeros(dt_number, 1);
dy = zeros(dt_number, 1);
fa0 = zeros(dt_number, 1);

n = 1;
for nn = 1:n_sensors - 1
   for mm = nn + 1:n_sensors
       dt(n) = base_points(nn, 4) - base_points(mm, 4);
       dr(n) = dt(n) * velocity;
       dx(n) = base_points(nn, 1) - base_points(mm, 1);
       dy(n) = base_points(nn, 2) - base_points(mm, 2);
       n = n + 1;
   end    
end

Y = -a_rad;
Y_prev = a_rad;

do = 1;

while (do == 1)  
    
   for n = 1:dt_number
      dfda(n) = -dx(n) * sin(a_rad) + dy(n) * cos(a_rad);
      fa0(n) = dx(n) * cos(a_rad) + dy(n) * sin(a_rad) + dr(n);     
   end

   Y = (1 / (dfda' * dfda)) * dfda' * fa0;   
   a_rad = a_rad - Y;
   
   it_cnt = it_cnt + 1;
   
   if ((it_cnt < max_iterations) && (abs(Y - Y_prev) > precision_threshold))
       do = 1;
   else 
       do = 0;
   end
   
   Y_prev = Y;   
end

a_rad = Nav_wrap_2pi(a_rad);

end