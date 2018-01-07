function [ CDFig, PPFig, PSFig, FSFig, focused ] = postProcessing( simVars, newFocuseRad, useNewRthresh, ...
                                                                   plotChargeDistribution, plotProbParticles, ... 
                                                                   plotPhaseSpace, recordPhaseSpace, plotFullSim, plotNotFocused )
%inputs:
% simVars should contain:
% 1. trajectory results: X,Y,U,W
% 2. Device results: V, MSE, zGrid, rGrid, Rq, Zq, Rb, Zb, Ez, Er, Mbleft, Mbright, Q 
% 3. Simulation corresponding input variables
%output:
% CDFIG - charge distribution
% PPFIG - problematic Particles (vector of figs)
% PSFIG - phase-Space
% FSFIG - Full Sim 

% TODO: define the input arguments so all functions below will have correct
% inputs

save('simVars-tmp.mat', '-struct', 'simVars');
load('simVars-tmp.mat');
delete('simVars-tmp.mat');
clear params;
CDFig = 0; PPFig = 0; PSFig = 0; FSFig = 0;
%--------------------------%
%Lens Post Processing
%--------------------------%

deviceLength = repetitions*distanceBetweenElectrodes + sideOffset;

if (plotChargeDistribution)
    leftElectrodeLength  = deviceRadius-leftElectrodeRadius;
    rightElectrodeLength = deviceRadius-rightElectrodeRadius;
    leftChargesWeight    = leftElectrodeLength/(leftElectrodeLength+rightElectrodeLength);
    NLeft                = round(N*leftChargesWeight);
    CDFig = displayChargeDistribution( repetitions, Rq, Q, N, NLeft, VaLeft, VaRight );
end

%-----------------------------%
% trajectories Post Processing
%-----------------------------%
if ( useNewRthresh)
    exitRthresh = newFocuseRad/1e3;
end
[startBetaX,startBetaZ,startBetaY, startBeta, startGamma] = getBG(Vz, Vx, Vy, ones(1,numOfParticles));
lastZ = zGrid(1,end);

[trajLen, hit, focused, hitMask, focusedMask, notFocusedMask ] = analyseTraj( Z, X, Zq, deviceRadius, leftElectrodeRadius,...
                                        rightElectrodeRadius, electrodeProximityThresh, exitRthresh, lastZ, numOfParticles);
[endBetaX, endBetaZ, endBeta, endGamma] = getBG(Vz,Vx,trajLen);
 row = 1:size(Y,1); 
 idx = sub2ind(size(Y), row, trajLen);
 exitR = Y(idx);
 
 lowAxialVel   = min(Vz(:,1)); 
 highAxialVel  = max(Vz(:,1));
 lowRadialVel  = min(Vx(:,1));
 highRadialVel = max(Vx(:,1));
 
 entryRvec = Y(:,1);
 
 
if(plotProbParticles)
    displayProbParticle(Z, X, Vz, zGrid, rGrid, Ez, Er, leftElectrodeRadius, rightElectrodeRadius,q, deviceRadius, VaLeft, VaRight,...
                                      m, multiParticle, trajectoryLen, axialEntryVel, radialEntryVel, Rr, numOfParticles, Zq);
                             %empty are variable for not phase space simulation
end

if(plotPhaseSpace)
    PSFig = displayPhaseSpace( focusedMask, numOfParticles, exitR, Rr, startGamma, startBetaX, endGamma,...
                       endBetaX, notFocusedMask, true, params);
end

if(recordPhaseSpace)
    phaseSpaceVideo( Z(passRows, :), X(passRows,:), Vz(passRows,:), Vx(passRows,:), Vy(passRows,:), Rr(passRows),...
                     zGrid, startBetaX(passRows), startGamma(passRows), params);
end

%-----------------------------%
% Full Simulation
%-----------------------------%

if(plotFullSim)
    trajParams = creatTrajParamsString(q, m, numOfParticles, hitMask, focused, endBetaZ,...
                                    endBetaX, exitR, electrodeProximityThresh, exitRthresh, lowAxialVel,...
                                    highAxialVel, lowRadialVel,highRadialVel, entryRvec, false, 0); 
                                
    fullParams = creatFullParamsString( electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                               deviceRadius, distanceBetweenElectrodes, deviceLength, ...
                                               VaLeft, VaRight, MSE, M, N, true, trajParams );
                                           
    FSFig = displayFullSim(VaLeft, VaRight, zGrid, rGrid, true, numOfParticles, ...
                                trajLen, V, X, Y, true, Mbleft, M, repetitions, ...
                                Zb, Rb, fullParams);
end


end

