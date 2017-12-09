function [ CDFig, PPFig, PSFig, FSFig, focused ] = postProcessing( simVars, newFocuseRad, useNewRthresh, ...
                                                                   plotChargeDistribution, plotProbParticles, ... 
                                                                   plotPhaseSpace, filmPhaseSpace, plotFullSim )
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
    exitRthresh = newFocuseRad;
end

[startBetaR,startBetaZ,startBeta, startGamma] = getBG(U, W, ones(1,numOfParticles));
lastZ = zGrid(1,end);
[trajLen, hit, focused, hitMask, focusedMask ] = analyseTraj( X, Y, Zq, deviceRadius, leftElectrodeRadius,...
                                        rightElectrodeRadius, electrodeProximityThresh, exitRthresh, lastZ, numOfParticles);
[endBetaR,endBetaZ,endBeta, endGamma] = getBG(U,W,trajLen);
 exitR = Y(trajLen);
 
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
                             trajectoryLen, axialEntryVelVec, radialEntryVelVec, entryRvec,...
                             [], [], [], numOfParticles); 
                             %empty are variable for not phase space simulation
end

if(plotPhaseSpace)
    PSFig = plotPhaseSpace( focusedMask, notFocusedMask, numOfParticles, exitR, entryRvec, ...
                            startGamma, startBetaR, endGamma, endBetaR, Y);
end

if(filmPhaseSpace)
    phaseSpaceVideo( q, m, numOfParticles, focused, focusedMask, X, Y, U, W, lowAxialVel, highAxialVel,...
                                lowRadialVel, highRadialVel, entryRvec, zGrid)
end

%-----------------------------%
% Full Simulation
%-----------------------------%

if(plotFullSim)
    trajParams = creatTrajParamsString(q, m, numOfParticles, hitMask, focused, endBetaZ,...
                                    endBetaR, exitR, proxTh, exitRthresh, lowAxialVel,...
                                    highAxialVel, lowRadialVel,highRadialVel, entryR); 
                                
    fullParams = creatParamsString( electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                               deviceRadius, distanceBetweenElectrodes, deviceLength, ...
                                               VaLeft, VaRight, MSE, M, N, simTraj, trajParams );
                                           
    FSFig = displayFullSim(VaLeft, VaRight, zGrid, rGrid, true, numOfParticles, ...
                                trajectoryLen, V, X, Y, true, Mbleft, M, repetitions, ...
                                Zb, Rb, fullParams);
end


end

