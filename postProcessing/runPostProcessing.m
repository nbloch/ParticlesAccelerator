function [ output_args ] = runPostProcessing( params )
%RUNPOSTPROCESSING Summary of this function goes here
%   Detailed explanation goes here

close all;
addpath(genpath([pwd]));
c0 = 3e8;
%load the simulation default values
M=730;
N=290;

load(sprintf('%s/simDefaultVals.mat',params.simPath));
save('defVal-tmp.mat', '-struct', 'defaultVals');
load('defVal-tmp.mat');
delete('defVal-tmp.mat');
clear defaultVals;

if(false)
elseif (params.distanceBetweenElectrodesSim)
    pdName = 'distanceBetweenElectrodes';
    pdVal = params.distanceBetweenElectrodesVal/(1e3);
    distanceBetweenElectrodes = params.distanceBetweenElectrodesVal/(1e3);
elseif (params.repetitionsSim)
    pdName = 'repetitions';
    pdVal = params.repetitionsVal;
    repetitions = params.repetitionsVal;
elseif (params.voltageSim)
    pdName = 'globalVa';
    pdVal = params.voltageVal*1e3;
    VaLeft = -params.voltageVal*1e3;
    VaRight = params.voltageVal*1e3;
elseif (params.electrodeRadiusSim)
    pdName = 'globalElectrodeRadius';
    pdVal = params.electrodeRadiusVal/(1e3);
    leftElectrodeRadius = params.electrodeRadiusVal/(1e3);
    rightElectrodeRadius = params.electrodeRadiusVal/(1e3);
else
    return
end

%load the simulation device's results
load(sprintf('%s/DeviceResults/%s-[%d].mat',params.simPath, pdName, pdVal));
%load the trajectories results
load(sprintf('%s/%s-[%d]-entryVel-[%d]/ParticleTrajectory.mat', params.simPath, pdName, pdVal, params.beta*c0));

parsStr=   {'V','MSE',...
         'zGrid', 'rGrid', 'Rq', 'Zq', 'Rb', 'Zb', 'Mbleft',...
         'Mbright', 'Ez', 'Er', 'Q','X','Y','U','W',...
         'beamInitialRadius', 'deviceRadius','distanceBetweenElectrodes',...
         'electrodeProximityThresh', 'electrodeWidth', 'exitRthresh',...
         'leftElectrodeRadius','numOfParticles', 'repetitions',...
         'rightElectrodeRadius', 'sideOffset', 'VaLeft' ,'VaRight', 'M', 'N', 'q', 'm'};

for i=1:length(parsStr(:))
    varname = parsStr{i};
    simVars.(varname)=eval(varname);
end



[ CDFig, PPFig, PSFig, FSFig, focused ] = postProcessing( simVars, params.newFocusRadius,  params.useNewFocusRadius, ...
                                                                   params.simulateChargeDist, params.simulateProbPart , ... 
                                                                   params.simulatePhaseSpace, params.simulatePhaseSpaceVideo,...
                                                                   params.simulateFullSim, params.plotNotFocused);

                                                               

end

