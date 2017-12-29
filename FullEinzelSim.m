function [ V, Z, X, Vz, Vx, MSE, fig, zGrid, rGrid, Rq, Zq, Rb, Zb, Ez, Er,...
    Mbleft, Mbright, NLeft, focused_particles, entryEmittance, exitEmittance, Q] =...
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
if (~plotChargeDistribution && ~plotPhaseSpaceVideo && ~plotProblematicParticles && ...
    ~plotFullSim            && ~plotEmittanceVsZ    && ~plotPhaseSpace           && ~dimensionsOptimization)
    fig = 0;
end
addpath(genpath(pwd));
%--------------------------------
%Params For Simulation
%--------------------------------
deviceLength = repetitions*distanceBetweenElectrodes + lensPreOffset + lensPostOffset;
q0 =  -1.60217662e-19;
eM = 9.10938356e-31;
if(simulateElectron)
    q = q0;
    m= eM;
else
    q = q*0;
    m = m*eM;
end
trajParams = {};
trajectoryLen = 0;
focused_particles = 0;
Z = []; X = []; Vz = []; Vx = [];

%--------------------------------
%Calculating Operating Point
%--------------------------------
if(dimensionsOptimization)
    tic;
    [ N, M, ~, fig.MSECalc ] = FindOperatingPoint(VaLeft, VaRight, electrodeWidth,leftElectrodeRadius,...
                                         rightElectrodeRadius, deviceRadius, distanceBetweenElectrodes,...
                                         use_bessel, repetitions, deviceLength,  convergeTh, growthTh );
    toc;
end

%--------------------------------
%Calculating Potential
%--------------------------------
if(calcPotential)
    tic;
    [V, zGrid, rGrid, MSE, Rq, Zq,  Rb, Zb, Mbleft, Mbright, NLeft, Q ] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                              deviceRadius, distanceBetweenElectrodes, N, M, repetitions, rPts, zPts, lensPreOffset, lensPostOffset);
    toc;
end
    
if (plotChargeDistribution); fig.ChargeDist = displayChargeDistribution( repetitions, Rq, Q, N, NLeft, VaLeft, VaRight ); end
deviceParams  = creatFullParamsString( electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                               deviceRadius, distanceBetweenElectrodes, deviceLength, ...
                                               VaLeft, VaRight, MSE, M, N);
%--------Process Results--------%
Rq = [Rq, -Rq];
Zq = [Zq,  Zq];
z_step = zGrid(1,2) - zGrid(1,1);
r_step = rGrid(2,1) - rGrid(1,1);
[Ez, Er] = gradient(-V, z_step, r_step);                                           
%--------------------------------
%Calculating Particle Trajectory
%--------------------------------
if (simulateMultiParticles)
    tic;
    [ Z, X, Vz, Vx, fig.PhaseSpace, fig.EmittanceVsZ, trajParams, trajectoryLen,...
                focused_particles, entryEmittance, exitEmittance] = ...
                 calcMultiParticles( V, zGrid, rGrid, q, m, Zq , electrodeProximityThresh,...
                 deviceRadius, leftElectrodeRadius, rightElectrodeRadius,...
                 numOfParticles, beamInitialRadius, maxInitMoment, Ek, exitRthresh,...
                 plotPhaseSpaceVideo, plotProblematicParticles, plotEmittanceVsZ,...
                 plotPhaseSpace, deviceParams);
    toc;
elseif (simulateSingleParticle)
    [  Z, X, Vz, Vx, trajParams ] = calcSingleParticle( V, zGrid, rGrid, q, m, axialEntryVelocity, radialEntryVelocity, entryR, exitRthresh,...
                                             Zq, Ez, Er, VaLeft, VaRight,...
                                             electrodeProximityThresh, deviceRadius, leftElectrodeRadius, rightElectrodeRadius, plotProblematicParticles);
end

%--------------------------------
%Plotting
%--------------------------------
if(plotFullSim)
    if(simulateMultiParticles || simulateSingleParticle); FullParams = [deviceParams ; trajParams];
    else; FullParams = deviceParams; end
    fig.FullSim = displayFullSim(VaLeft, VaRight, zGrid, rGrid, simulateMultiParticles, numOfParticles, ...
                                    trajectoryLen, V, Z, X, (simulateMultiParticles | simulateSingleParticle), Mbleft, M, repetitions, ...
                                    Zb, Rb, FullParams);
end         

end
