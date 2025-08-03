function [pos,vel] = simpsons_lunar_ephemeris(jd)

% Constants
tfac = 36525*3600*24;                  % seconds in a Julian century
t = (jd - 2451545.0) / 36525;           % centuries since J2000
% Amplitude matrix a (in km)
a = [383.0  31.5  10.6  6.2  3.2  2.3  0.8;
     351.0  28.9  13.7  9.7  5.7  2.9  2.1;
     153.2  31.5  12.5  4.2  2.5  3.0  1.8] * 1e3;
% Frequency matrix b (in rad/century)
b = [ 8399.685  70.990  16728.377  1185.622  7143.070  15613.745  8467.263;
      8399.687  70.997  8433.466  16728.380  1185.667  7143.058  15613.755;
      8399.672  8433.464 70.996  16728.364  1185.645  104.881  8399.116];
% Phase angle matrix c (in radians)
c = [5.381 6.169 1.453 0.481 5.017 0.857 1.010;
     3.811 4.596 4.766 6.165 5.164 0.300 5.565;
     3.807 1.629 4.595 6.162 5.167 2.555 6.248];
% Initialize vectors
pos = zeros(3, 1);
vel = zeros(3, 1);
% Compute position and velocity vectors
for i = 1:3
    for j = 1:7
        angle = b(i, j) * t + c(i, j);
        pos(i) = pos(i) + a(i, j) * sin(angle);
        vel(i) = vel(i) + a(i, j) * cos(angle) * b(i, j);
    end
    vel(i) = vel(i) / tfac;
end
% Display results
%fprintf('Moon Position Vector (km):\n');
%fprintf('Moon position vector at arrival (rm0) = [%.2f, %.2f, %.2f]\n', rm0(1), rm0(2), rm0(3));
%fprintf('moon distance at arrival: %.3f km\n', rm_mag);
%fprintf('Moon velocity vector at arrival (vm0) = [%.4f, %.4f, %.4f]\n', vm0(1), vm0(2), vm0(3));
%ym0 = [rm0(1), rm0(2), rm0(3), vm0(1), vm0(2), vm0(3)];
%fprintf('moon state vector at arrival in ECI frame (y0m) = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', ym0);
end
