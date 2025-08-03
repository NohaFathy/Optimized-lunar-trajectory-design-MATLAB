function [total_deltaV, total_fuel] = compute_total_deltavandfuel(alt0, mu_e, fac, vesc, Re, Isp, m0, g0)
    r0 = Re + alt0;
    v_circular = sqrt(mu_e / r0);
    v_TLI = fac * vesc;
    total_deltaV = v_TLI - v_circular;
    total_fuel = compute_fuel_used(total_deltaV, Isp, g0, m0);
end
