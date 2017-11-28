function [MSE] = getMSERepetitive(Q, M, N, VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                        deviceRadius, distanceBetweenElectrodes, use_bessel, repetitions )

    superM = 3*M;
    [r_q, z_q, r_b, z_b, Vm, ~, ~ ,~] = buildLens( superM, N, distanceBetweenElectrodes, electrodeWidth, deviceRadius,...
                                leftElectrodeRadius, rightElectrodeRadius, VaLeft, VaRight, repetitions);

    if(use_bessel)
        F = boundariesBesselPotential(r_b, z_b , r_q, z_q);
    else
        F = analytic_ring_green(r_b, z_b , r_q, z_q);
    end

    MSE = sum((F*Q - Vm).^2)/sum(Vm.^2);

end
