function [ Z, X, Vz, Vx, Vy, phaseSpaceFig, numErrFig, emittanceVsZFig, trajStr, ...
                 trajectoryLen, focused_particles, entryEmittance, exitEmittance, numericErr] = ...
                 calcMultiParticles (Voltage, zGrid, rGrid, q, m, Zq , proxTh, deviceRadius,...
                 leftElectrodeRadius, rightElectrodeRadius,...
                 numOfParticles, initialRadius, maxInitMoment, Ek, exitRthresh,...
                 recordPhaseSpace, plotProblematicParticles, plotEmittanceVsZ,...
                 plotPhaseSpace, deviceParams)
%------------------------------------%
%-----------Multi Particles----------%
%------------------------------------%                                         

%--------Constants---------%
c0 = 3e8;
electordesZ = unique(Zq);
entryZ   = zGrid(1,1);

%--------Init Variables---------%
phaseSpaceFig = 0;
emittanceVsZFig = 0;

focused = 0;
% particlesFailed = 0;

endBetaX = zeros(1,numOfParticles);
endBetaZ = zeros(1,numOfParticles);
% endBeta  = zeros(1,numOfParticles);
endGamma = zeros(1,numOfParticles);
exitR    = zeros(1,numOfParticles);

trajMaxLen = 100000;

Z =  zeros(numOfParticles,trajMaxLen);
X =  zeros(numOfParticles,trajMaxLen);
Vz = zeros(numOfParticles,trajMaxLen);
Vy = zeros(numOfParticles,trajMaxLen);
Vx = zeros(numOfParticles,trajMaxLen);
trajectoryLen = zeros(1,numOfParticles);

hitMask = zeros(1, numOfParticles);
focusedMask  = zeros(1,numOfParticles);
notFocusedMask = zeros(1,numOfParticles);

%--------Randomize Start Conditions---------%
%Vector of the angles in the XY plane (and not in relative to Z axis)
Rr = initialRadius*rand(1,numOfParticles);
phiR = pi*(-1+2*rand(1,numOfParticles));
Rx = Rr.*cos(phiR);
Ry = Rr.*sin(phiR);
Rr = Rr.*sign(Rx);

pr = maxInitMoment*rand(1, numOfParticles);
phiM = pi*(-1+2*rand(1, numOfParticles));
rPx = pr.*cos(phiM);
rPx(rPx.*Rx>0) = -rPx(rPx.*Rx>0);
rPy = pr.*sin(phiM);
rPy(rPy.*Ry>0) = -rPy(rPy.*Ry>0);

gamma0 = 1 + ((Ek*abs(q))/(m*(c0^2)));
rGamma = gamma0 *(0.999+0.002*rand(1,numOfParticles));

rPz = sqrt((rGamma.^2)-1-(rPx.^2)-(rPy.^2));

% figure()
% plot(Rr, rPx, 'o')


startBetaR = rPx./rGamma;
startBetaY = rPy./rGamma;
startGamma = rGamma;
Vyi = (rPy*c0)./rGamma;
%------------------------------------%
%--------Calculate Trajectory--------%
%------------------------------------%
parfor i = 1:numOfParticles
    [ trajectoriesZ,trajectoriesX, velocityZ, velocityX, hitElectrode ] =  relativeParticleTrajectory(Voltage, zGrid, rGrid, q, m, rPz(i)*c0/rGamma(i),...
                                                 rPx(i)*c0/rGamma(i), Rr(i), entryZ, proxTh,...
                                                 deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);
    if (hitElectrode)
        hitMask(i) = 1;
    else
        if ((abs(trajectoriesX(end)) <= exitRthresh) && (trajectoriesZ(end) >= zGrid(1,end)))
            focused = focused +1;
            focusedMask(i) = 1;    %exited the lens and beneath the focus threshold
        elseif ((trajectoriesZ(end) >= zGrid(1,end)) && (abs(trajectoriesX(end)) > exitRthresh))
            notFocusedMask(i) = 1; %exited the lens but exceeded the focus threshold
        end
        endBetaX(i) = velocityX(end)/c0;
        endBetaZ(i) = velocityZ(end)/c0;
        endGamma(i) = 1/sqrt(1-(endBetaX(i)^2+endBetaZ(i)^2+startBetaY(i)^2));
        exitR(i) = trajectoriesX(end);
    end
    len = length(trajectoriesZ);
    Z (i,:) = [trajectoriesZ,  NaN*ones(1,trajMaxLen-len)];
    X (i,:) = [trajectoriesX,  NaN*ones(1,trajMaxLen-len)];
    Vz (i,:) = [velocityZ,  NaN*ones(1,trajMaxLen-len)];
    Vx (i,:) = [velocityX,  NaN*ones(1,trajMaxLen-len)];
    Vy(i,:) = [Vyi(i)*ones(1,len), NaN*ones(1,trajMaxLen-len)]; 
    trajectoryLen (i) = length(trajectoriesZ); 
%         fprintf('Done Calculating Trajectory %d\n', i)
end

focusedMask = logical(focusedMask);
notFocusedMask = logical(notFocusedMask);
passMask = logical(focusedMask+notFocusedMask);
passRows = 1:size(Z,1);
passExitIdxs= sub2ind(size(Z), passRows(passMask), trajectoryLen(passMask));
exitPx = endBetaX(passMask).*endGamma(passMask);
exitPz = endBetaZ(passMask).*endGamma(passMask);
[entryEmittance, ~] = getEmittance( Rr(passMask), rPx(passMask), rPy(passMask), rPz(passMask), rGamma(passMask) );
[exitEmittance, ~] = getEmittance(X(passExitIdxs), exitPx, rPy(passMask), exitPz, endGamma(passMask)); 

%------------------------------------%
%----------Figures and Video---------%
%------------------------------------%
trajStr = creatTrajParamsString(q, m, numOfParticles, hitMask, focused, exitR, proxTh, exitRthresh, Ek, max(rPx), max(exitPx), Rr, ...
                                entryEmittance, exitEmittance);
params = [deviceParams; trajStr];
% --------Problematic Particle -------%
if (min(min(Vz)) < 0 && plotProblematicParticles)    
    displayProbParticle(Z, X, Vz, zGrid, rGrid, Ez, Er, leftElectrodeRadius, rightElectrodeRadius,q, deviceRadius, VaLeft, VaRight,...
                                      m, multiParticle, trajectoryLen, axialEntryVel, radialEntryVel, Rr, numOfParticles, Zq);
end
%---------Phase Space Figure---------%
if(plotPhaseSpace)
    phaseSpaceFig = displayPhaseSpace( passMask, numOfParticles, exitR, Rr, startGamma, startBetaR, endGamma,...
                       endBetaX, notFocusedMask, false, params);
end
%---------Emittance Figure---------%
[ emittance, zPtsCalc, numericErr] = calcEmiitanceVsZ( Z(passRows, :), X(passRows,:), Vz(passRows,:), Vx(passRows,:), Vy(passRows,:), zGrid);
if(plotEmittanceVsZ)
    [numErrFig, emittanceVsZFig] = plotEmittance( zPtsCalc, emittance, numericErr);
end
%---------Phase Space Video---------%
if(recordPhaseSpace)                         
    phaseSpaceVideo( Z(passRows, :), X(passRows,:), Vz(passRows,:), Vx(passRows,:), Vy(passRows,:), Rr(passRows),...
                   zGrid, startBetaR(passRows), startGamma(passRows), params)
end

focused_particles = sum(passMask);
fprintf('Done Calculating Trajectories\n')
end
