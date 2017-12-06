function [ CDFig, PPFig, PSFig, FSFig, focused ] = postProcessing( input_args )
%POSTPROCESSING Summary of this function goes here
%   Detailed explanation goes here
%output:
% CDFIG - charge distribution
% PPFIG - problematic Particles (vector of figs)
% PSFIG - phase-Space
% FSFIG - Full Sim 

% TODO: define the input arguments so all functions below will have correct
% inputs


%--------------------------%
%Lens Post Processing
%--------------------------%

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
[startBetaR,startBetaZ,startBeta, startGamma] = getBG(U,W,ones(1,numOfParticles));
[trajLen, hit, focused, hitMask, focusedMask ] = analyseTraj( X, Y, U, W, Zq, deviceRadius, leftElectrodeRadius,...
                                        rightElectrodeRadius, electrodeProximityThresh, exitRthresh);
[endBetaR,endBetaZ,endBeta, endGamma] = getBG(U,W,trajLen);
 exitR = Y(trajLen);
 
 trajParams = creatTrajParamsString(q, m, numOfParticles, hitMask, focused, endBetaZ,...
                                    endBetaR, exitR, proxTh, exitRthresh, lowAxialVel,...
                                    highAxialVel, lowRadialVel,highRadialVel, entryR); 

if(plotProbParticles)
    PPFig = plotProbParticle(X, Y, U, zGrid, rGrid, Ez, Er, leftElectrodeRadius, rightElectrodeRadius,...
                             q, deviceRadius, VaLeft, VaRight, m, simulatePhaseSpace,...
                             trajectoryLen, axialEntryVelVec, radialEntryVelVec, entryRvec,...
                             axialEntryVelocity, radialEntryVelocity, entryR, numOfParticles);
end

if(plotPhaseSpace)
    PSFig = plotPhaseSpace( focusedMask, notFocusedMask, numOfParticles, exitR, entryRvec, ...
                            startGamma, startBetaR, endGamma, endBetaR, Y);
end

if(recordPhaseSpace)
    phaseSpaceVideo( q, m, numOfParticles, focused, focusedMask, X, Y, U, W, lowAxialVel, highAxialVel,...
                                lowRadialVel, highRadialVel, entryRvec, zGrid)
end


if(plotFullSim)
    trajParams = creatTrajParamsString(q, m, numOfParticles, hitMask, focused, endBetaZ,...
                                    endBetaR, exitR, proxTh, exitRthresh, lowAxialVel,...
                                    highAxialVel, lowRadialVel,highRadialVel, entryR); 
                                
    fullParams = creatParamsString( electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                               deviceRadius, distanceBetweenElectrodes, deviceLength, ...
                                               VaLeft, VaRight, MSE, M, N, simTraj, trajParams );
                                           
    FSFig = displayFullSim(VaLeft, VaRight, zGrid, rGrid, simulatePhaseSpace, numOfParticles, ...
                                trajectoryLen, V, X, Y, simulateTrajectory, Mbleft, M, repetitions, ...
                                Zb, Rb, fullParams);
end


end

