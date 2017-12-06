function [ V, X, Y, U, W, MSE, fig, zGrid, rGrid, Rq, Zq, Rb, Zb, Ez, Er, Mbleft, Mbright, focused_particles ] =...
    FullEinzelSim(simPars)
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
c0 = 3e8;

if (~simulatePhaseSpace)
    numOfParticles = 1;
end

% if(useAngle)
%     %For one particle
%     axialEntryVelocity  = entryVel*cos(deg2rad(entryAngle));
%     radialEntryVelocity = entryVel*sin(deg2rad(entryAngle));
%     
% %     %For many particles
% %     lowAxialVel   = entryVel*cos(deg2rad(entryAngle));
% %     highAxialVel  = entryVel;
% %     lowRadialVel  = (1e-6)*c0;
% %     highRadialVel = entryVel*sin(deg2rad(entryAngle));
% end


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
      [V, zGrid, rGrid, MSE, fig(2), Rq, Zq,  Rb, Zb, Mbleft, Mbright ] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, ...
                                                    leftElectrodeRadius, rightElectrodeRadius, deviceRadius, ...
                                                    distanceBetweenElectrodes, N, M, use_bessel, repetitions, deviceLength,...
                                                    rPts, zPts);
      toc;
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
                                             useRelativeTrajectory, V, zGrid, rGrid, q, m, entryVel,...
                                             entryAngle, axialEntryVelocity, radialEntryVelocity, entryR, Zq ,...
                                             electrodeProximityThresh, beamInitialRadius,...
                                             lowAxialVel, highAxialVel, lowRadialVel, highRadialVel, numOfParticles,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, VaLeft, VaRight, Ez, Er, ...
                                             recordPhaseSpace, useAngle, beamInitialRadius);
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
