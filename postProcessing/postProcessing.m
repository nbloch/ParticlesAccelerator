function [ CDFig, PPFig, PSFig, FSFig, focused ] = postProcessing( simVars, newFocuseRad, useNewRthresh, ...
                                                                   plotChargeDistribution, plotProbParticles, ... 
                                                                   dispPhaseSpace, filmPhaseSpace, plotFullSim, plotNotFocused )
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

[startBetaR,startBetaZ,startBeta, startGamma] = getBG(U, W, ones(1,numOfParticles));
lastZ = zGrid(1,end);
[trajLen, hit, focused, hitMask, focusedMask, notFocusedMask ] = analyseTraj( X, Y, Zq, deviceRadius, leftElectrodeRadius,...
                                        rightElectrodeRadius, electrodeProximityThresh, exitRthresh, lastZ, numOfParticles);
[endBetaR, endBetaZ, endBeta, endGamma] = getBG(U,W,trajLen);
 row = 1:size(Y,1); 
 idx = sub2ind(size(Y), row, trajLen);
 exitR = Y(idx);
 
 lowAxialVel   = min(U(:,1)); 
 highAxialVel  = max(U(:,1));
 lowRadialVel  = min(W(:,1));
 highRadialVel = max(W(:,1));
 
 axialEntryVelVec  = U(:,1);
 radialEntryVelVec = W(:,1);
 entryRvec = Y(:,1);
 
 
if(plotProbParticles)
    PPFig = plotProbParticle(X, Y, U, zGrid, rGrid, Ez, Er, leftElectrodeRadius, rightElectrodeRadius,...
                             q, deviceRadius, VaLeft, VaRight, m, true,...
                             trajLen, axialEntryVelVec, radialEntryVelVec, entryRvec,...
                             [], [], [], numOfParticles,Zq); 
                             %empty are variable for not phase space simulation
end

if(dispPhaseSpace)
    PSFig = plotPhaseSpace( focusedMask, hitMask, numOfParticles, exitR, entryRvec, ...
                            startGamma, startBetaR, endGamma, endBetaR, Y, notFocusedMask, plotNotFocused);
end

if(filmPhaseSpace)
    phaseSpaceVideo( q, m, numOfParticles, focused, focusedMask, X, Y, U, W, lowAxialVel, highAxialVel,...
                     lowRadialVel, highRadialVel, entryRvec, zGrid, startBetaR, startGamma)
end

%-----------------------------%
% Full Simulation
%-----------------------------%

if(plotFullSim)
    trajParams = creatTrajParamsString(q, m, numOfParticles, hitMask, focused, endBetaZ,...
                                    endBetaR, exitR, electrodeProximityThresh, exitRthresh, lowAxialVel,...
                                    highAxialVel, lowRadialVel,highRadialVel, entryRvec, false, 0); 
                                
    fullParams = creatFullParamsString( electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                               deviceRadius, distanceBetweenElectrodes, deviceLength, ...
                                               VaLeft, VaRight, MSE, M, N, true, trajParams );
                                           
    FSFig = displayFullSim(VaLeft, VaRight, zGrid, rGrid, true, numOfParticles, ...
                                trajLen, V, X, Y, true, Mbleft, M, repetitions, ...
                                Zb, Rb, fullParams);
end


end

