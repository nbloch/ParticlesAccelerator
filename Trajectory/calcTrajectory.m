function [ Z, X, Vz, Vx, fig, parameters, trajectoryLen, focused_particles] = calcTrajectory( simulatePhaseSpace,...
                                              Voltage, zGrid, rGrid, q, m, entryVel,...
                                             entryAngle, axialEntryVelocity, radialEntryVelocity, entryR, Zq ,...
                                             electrodeProximityThresh, exitRthresh,...
                                             numOfParticles,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius,  ...
                                             recordPhaseSpace, useAngle, initialRadius, maxInitMoment, ...
                                             Ek)

c0 = 3e8;
e0 =  -1.60217662e-19;
eM = 9.10938356e-31;
electordesZ = unique(Zq);
entryZ   = zGrid(1,1);
trajectoryLen = 0;

if (simulatePhaseSpace)
    %------------------------------------%
    %-----------Multi Particles----------%
    %------------------------------------%
    focused = 0;
    particlesFailed = 0;

    endBetaR = zeros(1,numOfParticles);
    endBetaZ = zeros(1,numOfParticles);
    endBeta  = zeros(1,numOfParticles);
    endGamma = zeros(1,numOfParticles);
    exitR    = zeros(1,numOfParticles);

    trajMaxLen = 100000;
    
    Z =  zeros(numOfParticles,trajMaxLen);
    X =  zeros(numOfParticles,trajMaxLen);
    Vz = zeros(numOfParticles,trajMaxLen);
    Vx = zeros(numOfParticles,trajMaxLen);

    
    trajectoryLen = zeros(1,numOfParticles);

    entryRvec = initialRadius * rand(1, numOfParticles);
    %Vector of the angles in the XY plane (and not in relative to Z axis)
    phiVec = 2*pi*(rand(1, numOfParticles)-0.5);
    entryRvecSigned = entryRvec * sign(sin(phiVec));
    prVec = maxInitMoment*rand(1,numOfParticles);
    prVec(prVec.*entryRvecSigned >0) = -prVec(prVec.*entryRvecSigned >0);
    gamma0 = 1 + Ek/(eM*c0^2);
    rGamma = gamma0 *(0.999+rand(1,numOfParticles)/500);
    pzVec = sqrt(rGamma.^2 -1 -prVec.^2);
    
    
    startBetaR = prVec/(eM*c0);
    startBetaZ = pzVec/(eM*c0);
    startBeta  = sqrt(startBetaZ.^2+startBetaR.^2);
    startGamma = 1./sqrt(1-startBeta.^2);

    hitMask = zeros(1, numOfParticles);
    focusedMask  = zeros(1,numOfParticles);
    notFocusedMask = zeros(1,numOfParticles);
    %------------------------------------%
    %--------Calculate Trajectory--------%
    %------------------------------------%
    parfor i = 1:numOfParticles
        [ trajectoriesZ,trajectoriesX, velocityZ, velocityX, hitElectrode ] =  relativeParticleTrajectory(Voltage, zGrid, rGrid, q, m, pzVec(i)/eM,...
                                                     prVec(i)/eM, entryRvecSigned(i), entryZ, electrodeProximityThresh,...
                                                     deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);

        if (hitElectrode)
            particlesFailed = particlesFailed + 1;
            hitMask(i) = 1;
        else
            if ((abs(trajectoriesX(end)) <= exitRthresh) && (trajectoriesZ(end) >= zGrid(1,end)))
                focused = focused +1;
                focusedMask(i) = 1;
            elseif ((trajectoriesZ(end) >= zGrid(1,end)) && (abs(trajectoriesX(end)) > exitRthresh))
                notFocusedMask(i) = 1;
            end
            endBetaR(i) = velocityX(end)/c0;
            endBetaZ(i) = velocityZ(end)/c0;
            endBeta(i)  = sqrt(endBetaR(i)^2+endBetaZ(i)^2);
            endGamma(i) = 1/sqrt(1-endBeta(i)^2);
            exitR(i) = trajectoriesX(end);
        end
        len = length(trajectoriesZ);
        Z (i,:) = [trajectoriesZ,  NaN*ones(1,trajMaxLen-len)];
        X (i,:) = [trajectoriesX,  NaN*ones(1,trajMaxLen-len)];
        Vz (i,:) = [velocityZ,  NaN*ones(1,trajMaxLen-len)];
        Vx (i,:) = [velocityX,  NaN*ones(1,trajMaxLen-len)];

        trajectoryLen (i) = length(trajectoriesZ); 
%         fprintf('Done Calculating Trajectory %d\n', i)
    end
    focusedMask = logical(focusedMask);
    notFocusedMask = logical(notFocusedMask);
    %------------------------------------%
    %----------Figures and Video---------%
    %------------------------------------%
    
    %--------Problematic Particle -------%
%     if (min(min(U)) < 0)    
%        plotProbParticle(X, Y, U, zGrid, rGrid, Ez, Er, leftElectrodeRadius, rightElectrodeRadius,...
%                         q, deviceRadius, VaLeft, VaRight, m, simulatePhaseSpace,...
%                         trajectoryLen, axialEntryVelVec, radialEntryVelVec, entryRvec,...
%                         axialEntryVelocity, radialEntryVelocity, entryR, numOfParticles, Zq);
%     end
   
    %---------Phase Space Figure---------%
    fig = plotPhaseSpace( focusedMask, hitMask, numOfParticles, exitR, entryRvecSigned, ...
                                    startGamma, startBetaR, endGamma, endBetaR, X ,notFocusedMask, true);
                                
    parameters = creatTrajParamsString(q, m, numOfParticles, hitMask, focused, endBetaZ,...
                                       endBetaR, exitR, electrodeProximityThresh, exitRthresh, lowAxialVel,...
                                       highAxialVel, lowRadialVel,highRadialVel, entryRvecSigned, useAngle, entryVel);
       %---------Phase Space Video---------%
    if(recordPhaseSpace) 
       phaseSpaceVideo( q, m, numOfParticles, focused, focusedMask, Z, X, Vz, Vx, lowAxialVel, highAxialVel,...
                                lowRadialVel, highRadialVel, entryRvecSigned, zGrid, startBetaR, startGamma);
    end
else
    %------------------------------------%
    %-----------Single Particle----------%
    %------------------------------------%
    tic
    
    if(useAngle)
    %For one particle
        axialEntryVelocity  = entryVel*cos(deg2rad(entryAngle));
        radialEntryVelocity = entryVel*sin(deg2rad(entryAngle));
    end
    
    [ Z, X, Vz, Vx, hitElectrode] =  relativeParticleTrajectory(Voltage, zGrid, rGrid, q, m, axialEntryVelocity,...
                                             radialEntryVelocity, entryR, entryZ, electrodeProximityThresh,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);

    toc
    %------------------------------------%
    %---------------Figures--------------%
    %------------------------------------%
    fig = 0; %tmp
    if (min(Vz) < 0) %if particle started a motion backwards          
        plotProbParticle;
    end
    focused = 0;
    focusedMask = 0;
    if (X(end)*1e3 < exitRthresh) && (hitElectrode == 0)
        focused = 1;
        focusedMask =1;
    end
    parameters = {' ';
           '-------Particle Parameters:------';
          ['q: ', num2str(q/e0),'[e_0]'];
          ['M: ', num2str(m/eM),'[e_M]'];
          ' ';
           '----Particle Entry Parameters----';
          ['V_z_-_i_n: ', num2str(axialEntryVelocity/c0),'[c]'];
          ['V_r_-_i_n: ', num2str(radialEntryVelocity/c0),'[c]'];
          ['R_i_n: ',num2str(entryR*1e3),'[mm]'];
          ' ';
           '----Particles Exit Parameters----';
          ['V_z_-_o_u_t: ', num2str(Z(end)/c0),'[c]'];
          ['V_r_-_o_u_t: ', num2str(Vx(end)/c0),'[c]'];
          ['R_o_u_t: ',num2str(X(end)*1e3),'[mm]'];
          ['Hit: ', num2str(hitElectrode)]
          ['Focused: ', num2str(focused)];
          };
end
focused_particles = sum(focusedMask);
fprintf('Done Calculating Trajectories\n')
end
