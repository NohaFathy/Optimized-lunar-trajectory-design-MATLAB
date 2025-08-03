function [r0, v0, y0,w0_hat] = launch_conditions2rv(Re, dec0, alt0, RAAN0, gamma0,v0_mag,rm0)
    % Calculate the r0 vector at Trans-Lunar Injection (TLI)
    % r0 is the initial position vector of lunar trajectory insertion
    % Convert angles from degrees to radians
    dec_rad = deg2rad(dec0);
    RAAN_rad = deg2rad(RAAN0);
    gamma_rad = deg2rad(gamma0);
    % Compute position vector
    r0_mag = alt0 + Re;
    r0_x = r0_mag * cos(dec_rad) * cos(RAAN_rad);
    r0_y = r0_mag * cos(dec_rad) * sin(RAAN_rad);
    r0_z = r0_mag * sin(dec_rad);
    r0 = [r0_x, r0_y, r0_z];  % row vector
    ur_hat = r0/norm(r0);
    w0_hat = cross(r0,rm0)/ norm(cross(r0,rm0));
    ut_hat = cross (w0_hat,ur_hat);
    v0 = v0_mag * sin(gamma_rad) * ur_hat + v0_mag * cos(gamma_rad) * ut_hat;
  
    % Print vector as [x, y, z] format
    %fprintf('initial position vector of Probe (r0) = [%.2f, %.2f, %.2f] km\n', r0);
    %fprintf('initial velocity vector of Probe (v0) = [%.4f, %.4f, %.4f] km/s\n', v0);
    y0 = [r0_x, r0_y, r0_z v0 ];
    %fprintf('initial state vector of Probe in ECI frame (y0_prb) = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', y0);
   
end
