function [Q, Rq_vec, Zq_vec, Rb_vec, Zb_vec, MSE, NLeft, NRight, lastZoff,  MLeft, MRight] = getChargesRepetitive(VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                              deviceRadius, distanceBetweenElectrodes, N, M, use_bessel, repetitions )

    [r_q, z_q, r_b, z_b, Vm, NLeft, NRight, lastZoff, MLeft, MRight] = buildLens( M, N, distanceBetweenElectrodes, electrodeWidth, deviceRadius,...
                    leftElectrodeRadius, rightElectrodeRadius, VaLeft, VaRight, repetitions);

    if(use_bessel)
        F = boundariesBesselPotential(r_b, z_b , r_q, z_q);
    else
        F = analytic_ring_green(r_b, z_b , r_q, z_q);
    end

    Rq_vec = r_q(1,:);
    Zq_vec = z_q(1,:);
    Rb_vec = r_b(:,1);
    Zb_vec = z_b(:,1);
    Q = ((F.'*F)\(F.'))*Vm;
    MSE = getMSERepetitive(Q, M, N, VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                        deviceRadius, distanceBetweenElectrodes, use_bessel, repetitions );

end
