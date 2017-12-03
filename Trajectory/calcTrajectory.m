function [ X, Y, U, W, fig, parameters, trajectoryLen, focused_particles] = calcTrajectory( simulatePhaseSpace,...
                                             useRelativeTrajectory, Voltage, zGrid, rGrid, q, m, entryVel,...
                                             entryAngle, axialEntryVelocity, radialEntryVelocity, entryR, Zq ,...
                                             electrodeProximityThresh, exitRthresh,...
                                             lowAxialVel, highAxialVel, lowRadialVel, highRadialVel, numOfParticles,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, VaLeft, VaRight, Ez, Er, ...
                                             recordPhaseSpace, useAngle)

% if(useAngle)    
% %     %For many particles
% %     lowAxialVel   = entryVel*cos(deg2rad(entryAngle));
% %     highAxialVel  = entryVel;
% %     lowRadialVel  = (1e-6)*c0;
% %     highRadialVel = entryVel*sin(deg2rad(entryAngle));
% end

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

    X =  zeros(numOfParticles,20000);
    Y =  zeros(numOfParticles,20000);
    U = zeros(numOfParticles,20000);
    W = zeros(numOfParticles,20000);
    
    trajectoryLen = zeros(1,numOfParticles);
    if(useAngle)
        %EXPLANATION:
        %alpha = atan(((+-)electrodeRadius- Rin)/sideOffset
        %d=sideOffset
        entryRvec = (rand(1, numOfParticles) - 0.5)*3*leftElectrodeRadius;
        d = abs(zGrid(1) - electordesZ(1));
%         entryAngles = zeros(1, numOfParticles);
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
            if ((abs(trajectoriesY(end)) < exitRthresh) && (trajectoriesX(end) >= zGrid(1,end)))
                focused = focused +1;
                focusedMask(i) = 1;
            end
            endBetaR(i) = velocityY(end)/c0;
            endBetaZ(i) = velocityX(end)/c0;
            endBeta(i)  = sqrt(endBetaR(i)^2+endBetaZ(i)^2);
            endGamma(i) = 1/sqrt(1-endBeta(i)^2);
            exitR(i) = trajectoriesY(end);
        end
        len = length(trajectoriesX);
        X (i,:) = [trajectoriesX,  NaN*ones(1,20000-len)];
        Y (i,:) = [trajectoriesY,  NaN*ones(1,20000-len)];
        U (i,:) = [velocityX,  NaN*ones(1,20000-len)];
        W (i,:) = [velocityY,  NaN*ones(1,20000-len)];
        trajectoryLen (i) = length(trajectoriesX); 
        fprintf('Done Calculating Trajectory %d\n', i)
    end
    focusedMask = logical(focusedMask);

    %------------------------------------%
    %----------Figures and Video---------%
    %------------------------------------%
    
    %--------Problematic Particle -------%
%     if (min(min(U)) < 0)    
%             plotProbParticle
%     end
%    
    %---------Phase Space Figure---------%
    passMask = ~hitMask;
    notFocusedMask = logical((passMask)-focusedMask);
    entryPhase = startGamma.*startBetaR;
    endPhase   = endGamma.*endBetaR;
    
    fig = figure();
    legPlot(1) = plot(NaN,NaN,'r+');    hold on;
    legPlot(2) = plot(NaN,NaN,'b+');
    legPlot(3) = plot(NaN,NaN, '-');
    legPlot(4) = plot(NaN,NaN,'--');
    
    plot(exitR(focusedMask),  endPhase(focusedMask), 'r+')
    plot(entryRvec(focusedMask),entryPhase(focusedMask), 'b+')
    plot(exitR(notFocusedMask), endPhase(notFocusedMask), 'r+')
    plot(entryRvec(notFocusedMask), entryPhase(notFocusedMask), 'b+')
    
 
    i=1;
    while i <= numOfParticles
            while ((focusedMask(i) == 0) && (notFocusedMask(i) == 0) && i<= numOfParticles)
                i=i+1;
                if(i>numOfParticles); break; end
            end
            if (i>numOfParticles); break; end
            if (focusedMask(i) == 1)
                plot([exitR(i), entryRvec(i)],[endPhase(i), entryPhase(i)])
            else
                plot([exitR(i), entryRvec(i)],[endPhase(i), entryPhase(i)],'--')
            end
            i=i+1;
    end
    hold off
    leg = legend (legPlot,'Exit Phase', 'Start Phase', 'Focused', 'Not Focused');
    leg.Location = 'northeast';
    xl = max(abs([min(min(min(Y(passMask,:))),     min(entryRvec (passMask)))  max(max(max(Y(passMask,:))),     max(entryRvec(passMask)))]));
    yl = max(abs([min(min(min(endPhase(passMask))), min(entryPhase(passMask)))  max(max(max(endPhase(passMask))), max(entryPhase(passMask)))]));
    xlim(1.1.*[-xl xl]);
    ylim(1.1.*[-yl yl]);     
    xlabel('Entry/Exit R')
    ylabel('\gamma\beta_r')
    title('Particle Phase Space')
    
    %-------Parameters For Full Sim-------%
    exitLowAxialVel   = min(endBetaZ(~hitMask));
    exitHighAxialVel  = max(endBetaZ(~hitMask));
    exitLowRadialVel  = min(abs(endBetaR(~hitMask)));
    exitHighRadialVel = max(abs(endBetaR(~hitMask)));
    exitRRange        = max(abs(exitR(~hitMask)));

    parameters = {' ';
       '-------Particle Parameters-------';
      ['q: ', num2str(q/e0),'[e_0]'];
      ['M: ', num2str(m/eM),'[e_M]'];
      ['Number of Particles: ', num2str(numOfParticles)];
      ' ';
       '----Particle Entry Parameters----';
      ['V_z_-_i_n: [', num2str(lowAxialVel/c0),', ',num2str(highAxialVel/c0),'][c]'];
      ['V_r_-_i_n: [', num2str(lowRadialVel/c0),', ',num2str(highRadialVel/c0),'][c]'];
      ['R_i_n: [', num2str(-entryR*1e3),', ',num2str(entryR*1e3),'][mm]'];
      ' ';
       '----Particles Exit Parameters----';
      ['V_z_-_o_u_t: [', num2str(exitLowAxialVel),', ',num2str(exitHighAxialVel),'][c]'];
      ['V_r_-_o_u_t: [', num2str(exitLowRadialVel),', ',num2str(exitHighRadialVel),'][c]'];
      ['R_o_u_t: [', num2str(-exitRRange*1e3),', ',num2str(exitRRange*1e3),'][mm]'];
      ['d_p_r_o_x: ', num2str(electrodeProximityThresh*1e6), '[\mum]']
      ['R_f_o_c_u_s: ', num2str(exitRthresh*1e3), '[mm]'];
      ['Hit: ', num2str(particlesFailed)];
      ['Focused: ', num2str(focused)];
      };
       %---------Phase Space Video---------%
        gamma = 1./sqrt(1-((W.^2+U.^2)./c0^2));
        phase = gamma.*W./c0;
        startPhase = startGamma.*startBetaR;
        focusedMask = logical(focusedMask);
      if(recordPhaseSpace) 
          phaseSpaceVideo
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
