function fuel = fuel_used_fac(fac, alt0, mu_e, vesc, Re, Isp, m0, g0)
    r0 = Re + alt0;
    v_circular = sqrt(mu_e / r0);
    v_TLI = fac * vesc;
    deltaV = v_TLI - v_circular;

    fuel = compute_fuel_used(deltaV, Isp, g0, m0);
end
