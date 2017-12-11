function [ X, Y, U, W, fig, parameters, trajectoryLen, focused_particles] = calcTrajectory( simulatePhaseSpace,...
                                             useRelativeTrajectory, Voltage, zGrid, rGrid, q, m, entryVel,...
                                             entryAngle, axialEntryVelocity, radialEntryVelocity, entryR, Zq ,...
                                             electrodeProximityThresh, exitRthresh,...
                                             lowAxialVel, highAxialVel, lowRadialVel, highRadialVel, numOfParticles,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, VaLeft, VaRight, Ez, Er, ...
                                             recordPhaseSpace, useAngle, beamInitialRadius)

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
    
    X =  zeros(numOfParticles,trajMaxLen);
    Y =  zeros(numOfParticles,trajMaxLen);
    U = zeros(numOfParticles,trajMaxLen);
    W = zeros(numOfParticles,trajMaxLen);

    
    trajectoryLen = zeros(1,numOfParticles);
    if(useAngle)
        %EXPLANATION:
        %alpha = atan(((+-)electrodeRadius- Rin)/sideOffset
        %d=sideOffset
        entryRvec = (rand(1, numOfParticles) - 0.5)*2*beamInitialRadius;
        d = abs(zGrid(1) - electordesZ(1));
        upAngleBound = atan((leftElectrodeRadius-entryRvec)/d);
        downAngleBound= atan((-leftElectrodeRadius-entryRvec)/d);
        entryAngles = upAngleBound + (downAngleBound-upAngleBound).* rand(1, numOfParticles);
        axialEntryVelVec  = entryVel*cos(entryAngles);
        radialEntryVelVec = entryVel*sin(entryAngles);
    else
        entryR = abs(entryR);
        entryRvec         = (2*rand(1, numOfParticles)-1)*entryR;
        axialEntryVelVec  = rand(1, numOfParticles)*(highAxialVel  - lowAxialVel)  + lowAxialVel;
        radialEntryVelVec = rand(1, ceil(numOfParticles/2))*(highRadialVel - lowRadialVel) + lowRadialVel;
        radialEntryVelVec = [radialEntryVelVec, -(rand(1, floor(numOfParticles/2))*(highRadialVel - lowRadialVel))-lowRadialVel];
    end 
    
    startBetaR = radialEntryVelVec/c0;
    startBetaZ = axialEntryVelVec/c0;
    startBeta  = sqrt(startBetaZ.^2+startBetaR.^2);
    startGamma = 1./sqrt(1-startBeta.^2);

    hitMask = zeros(1, numOfParticles);
    focusedMask  = zeros(1,numOfParticles);
    notFocusedMask = zeros(1,numOfParticles);
    %------------------------------------%
    %--------Calculate Trajectory--------%
    %------------------------------------%
    parfor i = 1:numOfParticles
        if (useRelativeTrajectory)
        [ trajectoriesX,trajectoriesY, velocityX, velocityY, hitElectrode ] =  relativeParticleTrajectory(Voltage, zGrid, rGrid, q, m, axialEntryVelVec(i),...
                                                     radialEntryVelVec(i), entryRvec(i), entryZ, electrodeProximityThresh,...
                                                     deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);
        else
        [ trajectoriesX, trajectoriesY, velocityX, velocityY, hitElectrode ] =  ParticleTrajectory(Voltage, zGrid, rGrid, q, m, axialEntryVelVec(i),...
                                             radialEntryVelVec(i), entryRvec(i), entryZ, electrodeProximityThresh,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);
        end
        if (hitElectrode)
            particlesFailed = particlesFailed + 1;
            hitMask(i) = 1;
        else
            if ((abs(trajectoriesY(end)) <= exitRthresh) && (trajectoriesX(end) >= zGrid(1,end)))
                focused = focused +1;
                focusedMask(i) = 1;
            elseif ((trajectoriesX(end) >= zGrid(1,end)) && (abs(trajectoriesY(end)) > exitRthresh))
                notFocusedMask(i) = 1;
            end
            endBetaR(i) = velocityY(end)/c0;
            endBetaZ(i) = velocityX(end)/c0;
            endBeta(i)  = sqrt(endBetaR(i)^2+endBetaZ(i)^2);
            endGamma(i) = 1/sqrt(1-endBeta(i)^2);
            exitR(i) = trajectoriesY(end);
        end
        len = length(trajectoriesX);
        X (i,:) = [trajectoriesX,  NaN*ones(1,trajMaxLen-len)];
        Y (i,:) = [trajectoriesY,  NaN*ones(1,trajMaxLen-len)];
        U (i,:) = [velocityX,  NaN*ones(1,trajMaxLen-len)];
        W (i,:) = [velocityY,  NaN*ones(1,trajMaxLen-len)];

        trajectoryLen (i) = length(trajectoriesX); 
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
    fig = plotPhaseSpace( focusedMask, hitMask, numOfParticles, exitR, entryRvec, ...
                                    startGamma, startBetaR, endGamma, endBetaR, Y ,notFocusedMask, true);
    
    parameters = creatTrajParamsString(q, m, numOfParticles, hitMask, focused, endBetaZ,...
                                       endBetaR, exitR, electrodeProximityThresh, exitRthresh, lowAxialVel,...
                                       highAxialVel, lowRadialVel,highRadialVel, entryRvec);
       %---------Phase Space Video---------%
    if(recordPhaseSpace) 
       phaseSpaceVideo( q, m, numOfParticles, focused, focusedMask, X, Y, U, W, lowAxialVel, highAxialVel,...
                                lowRadialVel, highRadialVel, entryRvec, zGrid, startBetaR, startGamma);
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
    
    if (useRelativeTrajectory)
    [ X, Y, U, W, hitElectrode] =  relativeParticleTrajectory(Voltage, zGrid, rGrid, q, m, axialEntryVelocity,...
                                             radialEntryVelocity, entryR, entryZ, electrodeProximityThresh,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);
    else
    [ X, Y, U, W ,hitElectrode] =  ParticleTrajectory(Voltage, zGrid, rGrid, q, m, axialEntryVelocity,...
                                             radialEntryVelocity, entryR, entryZ, electrodeProximityThresh,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);
    end
    toc
    %------------------------------------%
    %---------------Figures--------------%
    %------------------------------------%
    fig = 0; %tmp
    if (min(U) < 0) %if particle started a motion backwards          
        plotProbParticle;
    end
    focused = 0;
    focusedMask = 0;
    if (Y(end)*1e3 < exitRthresh) && (hitElectrode == 0)
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
          ['V_z_-_o_u_t: ', num2str(X(end)/c0),'[c]'];
          ['V_r_-_o_u_t: ', num2str(W(end)/c0),'[c]'];
          ['R_o_u_t: ',num2str(Y(end)*1e3),'[mm]'];
          ['Hit: ', num2str(hitElectrode)]
          ['Focused: ', num2str(focused)];
          };
end
focused_particles = sum(focusedMask);
fprintf('Done Calculating Trajectories\n')
end
