function [ V, X, Y, U, W, MSE, fig, zGrid, rGrid, Rq, Zq, Rb, Zb, Ez, Er, Mbleft, Mbright, focused_particles, Q ] =...
    FullEinzelSim(simPars) %#ok<INUSD>
%NOTE: The parameters "calcPotential, V, zGrid, rGrid, Rq, Zq, MSE" are
% used for device potential caching. if calcPotential is set to true, the
% other parameters are unused and can be left empty

%Load struct parameters into workspace
save('params-tmp.mat', '-struct', 'simPars');
load('params-tmp.mat');
delete('params-tmp.mat');
clear params;
%END

close all;
fig = gobjects(1,4);
addpath(genpath([pwd]));
%--------------------------------
%Params For trajectory
%--------------------------------

deviceLength = repetitions*distanceBetweenElectrodes + sideOffset;

if(simulateElectron)
   q =  -1.60217662e-19;
   m = 9.10938356e-31;
end

if (~simulatePhaseSpace)
    numOfParticles = 1;
end


%--------------------------------
%Calculating Operating Point
%--------------------------------
if(dimensionsOptimization)
    tic;
    [ N, M, ~, fig(1) ] = FindOperatingPoint(VaLeft, VaRight, electrodeWidth,leftElectrodeRadius,...
                                         rightElectrodeRadius, deviceRadius, distanceBetweenElectrodes,...
                                         use_bessel, repetitions, deviceLength,  convergeTh, growthTh );
    toc;
end


%--------------------------------
%Calculating Potential
%--------------------------------
if(calcPotential)
    tic;
    [V, zGrid, rGrid, MSE, fig(2), Rq, Zq,  Rb, Zb, Mbleft, Mbright, Q ] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, ...
                                                leftElectrodeRadius, rightElectrodeRadius, deviceRadius, ...
                                                distanceBetweenElectrodes, N, M, use_bessel, repetitions, deviceLength,...
                                                rPts, zPts);
    toc;
else
    Q=0;
end
%--------------------------------
%Calculating Particle Trajectory
%--------------------------------
Rq = [Rq, -Rq];
Zq = [Zq,  Zq];
X = []; Y = []; U = []; W = [];
z_step = zGrid(1,2) - zGrid(1,1);
r_step = rGrid(2,1) - rGrid(1,1);
[Ez, Er] = gradient(-V, z_step, r_step);
if(simulateTrajectory)
  tic;
  [ X, Y, U, W, fig(3),trajParams, trajectoryLen, focused_particles] = calcTrajectory( simulatePhaseSpace,...
                                              V, zGrid, rGrid, q, m, entryVel,...
                                             entryAngle, axialEntryVelocity, radialEntryVelocity, entryR, Zq ,...
                                             electrodeProximityThresh, beamInitialRadius,...
                                             numOfParticles,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, ...
                                             recordPhaseSpace, useAngle, beamInitialRadius, maxInitMoment, Ek);
  toc;
end

%--------------------------------
%Plotting
%--------------------------------
FullParams  = creatFullParamsString( electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                               deviceRadius, distanceBetweenElectrodes, deviceLength, ...
                                               VaLeft, VaRight, MSE, M, N, simulateTrajectory, trajParams );

fig(4) = displayFullSim(VaLeft, VaRight, zGrid, rGrid, simulatePhaseSpace, numOfParticles, ...
                                trajectoryLen, V, X, Y, simulateTrajectory, Mbleft, M, repetitions, ...
                                Zb, Rb, FullParams);
                            
end
