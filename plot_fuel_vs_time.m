function plot_fuel_vs_time(t, m0, deltaV_opt, burn_duration, Isp)
    g0 = 9.80665;  % m/s^2
    deltaV_mps = deltaV_opt * 1000;

    fuel_total = m0 * (1 - exp(-deltaV_mps / (Isp * g0)));

    fuel_profile = zeros(size(t));
    idx_burn = t <= burn_duration;

    
    fuel_profile(idx_burn) = linspace(0, fuel_total, sum(idx_burn));
    fuel_profile(~idx_burn) = fuel_total;

    figure;
    plot(t/3600, fuel_profile, 'r', 'LineWidth', 2);
    xlabel('Time (hours)');
    ylabel('Fuel Used (kg)');
    title('Main Maneuver Fuel Usage');
    grid on;
end
