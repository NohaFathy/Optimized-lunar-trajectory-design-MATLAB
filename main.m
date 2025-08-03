clear all; close all; clc
%% General data
m0 =  28801;  
g0 = 9.80665;  
Isp=   311;
deg = pi/180;
days = 24*3600;
Re= 6378;
Rm = 1737;
alpha0 =90;
gamma0 = 40;
m_e = 5974.e21;
m_m = 73.48e21;   
mu_e = 398600.4;
mu_m = 4902.8;
D = 384400;
RS = D*(m_m/m_e)^(2/5);
Title = 'Start';
fac = .9924;
ttt = 3*days;
tf = ttt+ 2.667*days;
dec0 = 15;
alt0 =320;
z0 = alt0+Re;
speed_at_LEO = sqrt(mu_e / z0); 
vesc = sqrt(2* mu_e / z0);
speed_TLI= fac * vesc;
delta_v = speed_TLI - speed_at_LEO;

% compute initial fuel
fuel_before_optimization = compute_fuel_used(delta_v, Isp, g0, m0);
efficiency = (fuel_before_optimization / m0) * 100;

% Output
fprintf('\n\n%s\n\n', Title)

fprintf('\n=== Before Optimization ===\n');
fprintf('fac (nominal)           = %.5f\n', fac);
fprintf('ΔV used (km/s)          = %.4f\n', delta_v);
fprintf('Fuel used (kg)          = %.4f\n', fuel_before_optimization);
fprintf('Efficiency (%%)         = %.2f %%\n', efficiency);

% === Optimization ===
fac_min = 0.97;
fac_max = 1.01;
objective = @(fac) fuel_used_fac(fac, alt0, mu_e, vesc, Re, Isp, m0, g0);
[fac_opt, min_fuel_used] = fminbnd(objective, fac_min, fac_max);
[deltaV_opt, ~] = compute_total_deltavandfuel(alt0, mu_e, fac_opt, vesc, Re, Isp, m0, g0);
opt_efficiency = (min_fuel_used / m0) * 100;

fprintf('\n=== Optimization Results ===\n');
fprintf('Optimal fac              = %.5f\n', fac_opt);
fprintf('ΔV (km/s)                = %.4f\n', deltaV_opt);
fprintf('Fuel used (kg)           = %.4f\n', min_fuel_used);
fprintf('Efficiency (%%)          = %.2f %%\n', opt_efficiency);

% === Plot fuel vs fac ===
fac_vals = linspace(0.97, 1.02, 200);
fuel_vals = arrayfun(@(f) compute_fuel_used(f*vesc - speed_at_LEO, Isp, g0, m0), fac_vals);

figure;
plot(fac_vals, fuel_vals, 'LineWidth', 2);
xlabel('fac (fraction of v_{esc})');
ylabel('Fuel used (kg)');
title('Fuel Used vs fac');
grid on;

%% Date and time of lunar arrival:
year = 2020;
month = 5;
day= 4;
hour = 12;
minute = 0;
second = 0;
UT =12;
t0 =0;
jd0 = julian_day(year, month, day, UT);

%% state vector of the moon at target date 

[rm0,vm0] = simpsons_lunar_ephemeris(jd0);
[RA, Dec] =  ra_and_dec_from_r(rm0);
distance =norm(rm0);
hmoon = cross(rm0,vm0);
hmoon_mag = norm(hmoon);
inclmoon = acosd(hmoon(3)/hmoon_mag);

fprintf('Date and time of arrival at moon: ')
fprintf('%s/%s/%s %s:%s:%s', ...
num2str(month), num2str(day), num2str(year), ...
num2str(hour), num2str(minute), num2str(second))
fprintf('\nMoon''s position: ')
fprintf(' distance             = %11g km\n', distance);
fprintf(' RA                   = %11g deg\n', RA);
fprintf(' Dec                  = %11g deg\n', Dec);
fprintf(' inclmoon             = %11g deg\n', inclmoon);

%% ...Initial position vector of probe:
I = [1;0;0];
J = [0;1;0];
K = cross(I,J);
r0_mag = Re+alt0;
r0 = r0_mag*(cosd(alpha0)*cosd(dec0)*I + ...
 sind(alpha0)*cosd(dec0)*J + ...
 sind(dec0)*K);

w0 =  cross(r0,rm0)/norm(cross(r0,rm0));

fprintf('\nThe probe at Earth departure (t = %g sec):\n', t0);
fprintf(' Altitude             = %11g km\n', alt0);
fprintf(' Right Ascension      = %11g deg\n', alpha0);
fprintf(' Declination          = %11g deg\n', dec0);
fprintf(' Flight Path Angle    = %11g deg\n', gamma0);
fprintf(' Speed                = %11g km/s\n', speed_TLI);
fprintf(' Escape Speed         = %11g km/s\n', vesc);
fprintf(' v / v_esc            = %11g\n', speed_TLI/vesc);
fprintf(' Inclination of Translunar Orbit = %11g deg\n', acosd(w0(3)));

 %% ...Initial velocity vector of probe:
 ur = r0/norm(r0);
 uperp = cross(w0,ur)/norm(cross(w0,ur));
 vr_mag =speed_TLI*sind(gamma0);
 vperp_mag = speed_TLI*cosd(gamma0);
 v0 = vr_mag*ur + vperp_mag*uperp;
 uv0 = v0 /speed_TLI;
%...Initial state vector of the probe:
 y0 = [r0(1) r0(2) r0(3) v0(1) v0(2) v0(3)]';

%% ODE integration
options = odeset('RelTol', 1.e-10, 'AbsTol', 1.e-10,'Stats', 'off');
rates_func = @(t, y) rates(t, y, jd0, ttt, days, mu_m, mu_e);
[t, y] = ode45(rates_func, [t0 tf], y0, options);
%[r0 ,v0, y0,w0_hat] = launch_conditions2rv(Re,dec0, alt0, RAAN0, gamma0,v0_mag,rm0);

%% ...Spacecraft trajectory
 % in ECI frame:
 X =y(:,1); Y = y(:,2); Z = y(:,3);
 vX = y(:,4); vY = y(:,5); vZ = y(:,6);
 % in Moon-fixed frame:
 x =[]; y =[]; z =[];
 %% ...Moon trajectory
 % in ECI frame:
 Xm = []; Ym = []; Zm = [];
 vXm = []; vYm = []; vZm = [];
 % in Moon-fixed frame:
 xm = []; ym = []; zm = [];

%Starting value in the search for perilune
dist_min = 1e30;

for i = 1:length(t)
    ti = t(i);

    %...Probe's inertial position vector at time ti
    r = [X(i); Y(i); Z(i)];

     %...Moon's inertial position and velocity vectors at time ti:
    jd = jd0 - (ttt - ti) / days;
    [rm_mag, vm] = simpsons_lunar_ephemeris(jd);

    %...Moon's inertial state vector at time ti:
    Xm = [Xm; rm_mag(1)];
    Ym = [Ym; rm_mag(2)];
    Zm = [Zm; rm_mag(3)];
    vXm = [vXm; vm(1)];
    vYm = [vYm; vm(2)];
    vZm = [vZm; vm(3)];

     %...Moon-fixed rotating xyz frame:
    x_axis = rm_mag;
    z_axis = cross(rm_mag, vm);
    y_axis = cross(z_axis, x_axis);
    i_hat = x_axis / norm(x_axis);
    j_hat = y_axis / norm(y_axis);
    k_hat = z_axis / norm(z_axis);

    %...DCM of transformation from ECI to moon-fixed frame:
    Q = [i_hat'; j_hat'; k_hat'];

    %...Components of probe's inertial position vector in moon-fixed frame:
    r_moon_fixed = Q * r;
    x = [x; r_moon_fixed(1)];
    y = [y; r_moon_fixed(2)];
    z = [z; r_moon_fixed(3)];

    %...Components of moon's inertial position vector in moon-fixed frame:
    rm_moon_fixed = Q * rm_mag;
    xm = [xm; rm_moon_fixed(1)];
    ym = [ym; rm_moon_fixed(2)];
    zm = [zm; rm_moon_fixed(3)];

     %...Find perilune of the probe:
    distance_vector = r - rm_mag;
    current_distance = norm(distance_vector);

    if current_distance < dist_min
        imin = i;
        dist_min = current_distance;
    end
end

 %...Location of the Moon at TLI:
rmTLI = [Xm(1); Ym(1); Zm(1)];
[RATLI, DecTLI] = ra_and_dec_from_r(rmTLI);
fprintf('\nThe moon when the probe is at TLI:')
fprintf(' Distance             = %11g km\n', norm(rmTLI));
fprintf(' Right Ascension      = %11g deg\n', RATLI);
fprintf(' Declination          = %11g deg\n', DecTLI);


 %...Spacecraft velocity at perilune:
v_atdmin =[vX(imin); vY(imin); vZ(imin) ];
 %...State vector and celestial position of moon when probe is at perilune:
rm_perilune = [Xm(imin); Ym(imin); Zm(imin)];
vm_perilune = [vXm(imin); vYm(imin); vZm(imin)];
[RA_at_perilune,Dec_at_perilune] = ra_and_dec_from_r(rm_perilune);
target_error = norm(rm_perilune - rm0);   

fprintf('\nThe moon when the probe is at perilune:\n');
fprintf(' Distance             = %11g km\n', norm(rm_perilune));
fprintf(' Speed                = %11g km/s\n', norm(vm_perilune));
fprintf(' Right Ascension      = %11g deg\n', RA_at_perilune);
fprintf(' Declination          = %11g deg\n', Dec_at_perilune);
fprintf(' Target Error         = %11g km\n', target_error);




 %...Speed of probe relative to Moon at perilune:
rel_speed = norm(v_atdmin - vm_perilune);

 %...End point of trajectory:
rend = [X(end); Y(end); Z(end)];
alt_end = norm(rend) - Re;
[ra_end, dec_end] = ra_and_dec_from_r(rend);

 %...Find the history of the trajectory's binormal:
for i = 1:imin
    time(i)     = t(i);
    
    r_vec       = [X(i); Y(i); Z(i)];
    r_mag       = norm(r_vec);
    
    v_vec       = [vX(i); vY(i); vZ(i)];
    
    rm_vec      = [Xm(i); Ym(i); Zm(i)];
    rm_mag      = norm(rm_vec);
    
    rms_vec     = rm_vec - r_vec;
    rms_mag     = norm(rms_vec);
    
    a_earth     = -mu_e * r_vec / r_mag^3;
    a_moon      = mu_m * (rms_vec / rms_mag^3 - rm_vec / rm_mag^3);
    
    atotal      = a_earth + a_moon;
    
    binormal    = cross(v_vec, atotal) / norm(cross(v_vec, atotal));
    binormal_z  = binormal(3);
    
    incl(i)     = acosd(binormal_z);
    moon_distance(i) = norm(r_vec - rm_vec);  % distance from the moon

end

figure;
plot(moon_distance, incl, 'b', 'LineWidth', 2);
hold on;

plot(moon_distance(1), incl(1), 'go', 'MarkerSize', 10, 'LineWidth', 2); % TLI
text(moon_distance(1), incl(1)+0.5, 'TLI', 'Color', 'green', 'FontSize', 12);

plot(moon_distance(imin), incl(imin), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Perilune
text(moon_distance(imin), incl(imin)+0.5, 'Perilune', 'Color', 'red', 'FontSize', 12);

xlabel('Distance from Moon (km)');
ylabel('Inclination (degrees)');
title('Variation of Inclination with Distance from the Moon');
grid on;
legend('Inclination profile','TLI','Perilune');


%%... Output
fprintf('\nThe probe at Perilune:\n');
fprintf(' Altitude above Moon        = %11.4f km\n', dist_min - Rm);
fprintf(' Absolute Speed             = %11.4f km/s\n', norm(v_atdmin));
fprintf(' Relative Speed to Moon     = %11.4f km/s\n', rel_speed);
fprintf(' Inclination of Osculating Plane = %11.4f deg\n', incl(imin));

fprintf('\nTime from TLI to Perilune = %11.4f hours (%.4f days)\n', abs(t(imin))/3600, abs(t(imin))/(3600*24));
fprintf('Total Time of Flight          = %11.4f days\n', t(end)/days);
fprintf('Time to Target Point          = %11.4f days\n', ttt/days);

fprintf('\nFinal Earth Altitude        = %11.4f km\n', alt_end);
fprintf('Final Right Ascension         = %11.4f deg\n', ra_end);
fprintf('Final Declination             = %11.4f deg\n', dec_end);
 %...End output


burn_duration = 1800;  
plot_fuel_vs_time(t, m0, deltaV_opt, burn_duration, Isp);



v = sqrt(vX.^2 + vY.^2 + vZ.^2);

figure;
plot(t/3600, v, 'b', 'LineWidth', 2);
xlabel('Time (hours)');
ylabel('Velocity (km/s)');
title('Velocity of Probe vs Time');
grid on;

%plotit_XYZ(X, Y, Z, Xm, Ym, Zm, imin, Re, Rm)
%plotitxyz(x, y, z, xm, ym, zm, imin,Re,Rm)

