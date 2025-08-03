function fuel = compute_fuel_used(deltaV, Isp, g0, m0)
    if deltaV <= 0
        fuel = 0;
        return;
    end
    deltaV_mps = deltaV * 1000;
    mf = m0 / exp(deltaV_mps / (Isp * g0));
    fuel = m0 - mf;
end
