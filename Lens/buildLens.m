function [r_q, z_q, r_b, z_b, Vm, NLeft, NRight, lastZoff, MLeft, MRight] = buildLens( M, N, distanceBetweenElectrodes, electrodeWidth, ...
                                                  deviceRadius, leftElectrodeRadius, rightElectrodeRadius,VaLeft, VaRight, repetitions)

                                               

    [ r_q, z_q , r_b, z_b, MLeft, MRight, NLeft, NRight] = createUnitCell( M, N, distanceBetweenElectrodes,...
                                    electrodeWidth, deviceRadius, leftElectrodeRadius, rightElectrodeRadius );

    Vm = [VaLeft*ones(MLeft,1); VaRight*ones(MRight,1)];
    if(repetitions >1)
        for i=1:repetitions-1
            [tmpRQ, tmpZQ, tmpRB, tmpZB, MLeft, MRight, NLeft, NRight] = createUnitCell( M, N, distanceBetweenElectrodes,...
                                    electrodeWidth, deviceRadius, leftElectrodeRadius, rightElectrodeRadius );

            tmpRQ = repmat(tmpRQ,i,1);
            tmpZQ = repmat(tmpZQ,i,1);
            tmpZB = repmat(tmpZB,1,i);
            tmpRB = repmat(tmpRB,1,i);
            tmpZQ = tmpZQ + 2*distanceBetweenElectrodes*i;
            tmpZB = tmpZB + 2*distanceBetweenElectrodes*i;

            r_q = [r_q, tmpRQ];
            r_q = [r_q; repmat(r_q(1,:),M,1)];
            z_q = [z_q, tmpZQ];
            z_q = [z_q; repmat(z_q(1,:),M,1)];

            r_b = [r_b; tmpRB];
            r_b = [r_b, repmat(r_b(:,1),1,N)];
            z_b = [z_b; tmpZB];
            z_b = [z_b, repmat(z_b(:,1),1,N)];

            Vm = [Vm; VaLeft*ones(MLeft,1); VaRight*ones(MRight,1)];

        end
    end

    [tmpRQ, tmpZQ, tmpRB, tmpZB, lastZoff] = createLastElectrode(MLeft, NLeft, distanceBetweenElectrodes, ...
                        electrodeWidth, deviceRadius, leftElectrodeRadius);

    tmpRQ = [tmpRQ; repmat(tmpRQ(end,:),MRight,1)];
    tmpZQ = [tmpZQ; repmat(tmpZQ(end,:),MRight,1)];
    tmpZB = [tmpZB, repmat(tmpZB(:,end),1,NRight)];
    tmpRB = [tmpRB, repmat(tmpRB(:,end),1,NRight)];

    tmpRQ = repmat(tmpRQ,repetitions,1);
    tmpZQ = repmat(tmpZQ,repetitions,1);
    tmpZB = repmat(tmpZB,1,repetitions);
    tmpRB = repmat(tmpRB,1,repetitions);
    tmpZQ = tmpZQ + 2*distanceBetweenElectrodes*(repetitions-1)+lastZoff;%-deviceLength;%+lastZoff;
    tmpZB = tmpZB + 2*distanceBetweenElectrodes*(repetitions-1)+lastZoff;%-deviceLength;%+lastZoff;

    r_q = [r_q, tmpRQ];
    r_q = [r_q; repmat(r_q(1,:),MLeft,1)];
    z_q = [z_q, tmpZQ];
    z_q = [z_q; repmat(z_q(1,:),MLeft,1)];

    r_b = [r_b; tmpRB];
    r_b = [r_b, repmat(r_b(:,1),1,NLeft)];
    z_b = [z_b; tmpZB];
    z_b = [z_b, repmat(z_b(:,1),1,NLeft)];

    %Centralization of the electrodes in the device
    centralizationFactor = (z_q(end, end) + z_q(1,1))/2;
    z_q = z_q-centralizationFactor;
    z_b = z_b-centralizationFactor;

    Vm = [Vm; VaLeft*ones(MLeft,1)];
end
