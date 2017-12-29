c0 = 3e8;
e0 = -1.60217662e-19;
eM = 9.10938356e-31;          

%device Parameters
params.deviceRadius              = 5/1e6;
params.lensPreOffset             = 2/(1e6);
params.lensPostOffset            = 2/(1e6);
params.distanceBetweenElectrodes = 1/1e6;
params.leftElectrodeRadius       = 1/1e6;
params.rightElectrodeRadius      = 1/1e6;
params.electrodeWidth            = 0.01/1e6;
params.repetitions               = 1;
params.VaLeft                    = -15;
params.VaRight                   = 15;

%Resolution Parameters     
params.M                         = 730;
params.N                         = 290;     
params.rPts                      = 500;
params.zPts                      = 500;
params.dimensionsOptimization    = false;
params.convergeTh                = 5e-8;
params.growthTh                  = 0.5e-4;

%what To Simulate
params.dontSimulateTrajectory  =  false;
params.simulateSingleParticle  =  false;
params.simulateMultiParticles  =  true;

%Figures
params.plotFullSim              = true;
params.plotPhaseSpace           = true;
params.plotPhaseSpaceVideo      = false;
params.plotChargeDistribution   = true;
params.plotProblematicParticles = true;
params.plotEmittanceVsZ         = true;

%Particle Parameters
params.electrodeProximityThresh  = 0.01/1e6;
params.q                         = 1;
params.m                         = 1;
params.simulateElectron          = true;

%Multi Particle
params.numOfParticles    = 1000;
params.beamInitialRadius = 0.5/1e6;
params.maxInitMoment     = 1e-3;
params.Ek                = 10*1e3;

%Single Particle - irrelevant for this sim
params.axialEntryVelocity  = 0.1*c0;
params.radialEntryVelocity = 0.05*c0;
params.entryR              = 0;
params.exitRthresh         = 0.25/1e6;

%Save Sim
params.simName       = 'lens';
params.simGlobalName = 'Simulations - Emittance Style6';
params.eraseOldSim   = true;
params.genNewSeed    = true;
params.savePlots     = true;
%params.simPdName     = 'Repetitions';
params.SingleSim     = false;

paramsFields = fieldnames(params);

iterParamsCharges = struct();
iterParamsDevice = struct();

%charge parameters iteration variables definitions

iterParamsCharges.Ek                        = [10e3, 100e3];    %eV units
iterParamsCharges.maxInitMoment             = [1e-3, 1e-2];


%device parameters iteration variables definitions

iterParamsDevice.globalVa = [10,15,30,50,100,300,500,750,1000,2000];
%iterParamsDevice.repetitions = [1,2,3,5,7,9];
iterParamsDevice.distanceBetweenElectrodes=[0.5,0.75,1,1.25,2]/1e6;
iterParamsDevice.globalElectrodeRadius=[0.25,0.5,0.75,1,1.5,2]/1e6;
iterParamsDevice.params.lensPreOffset = [1,1.5,2,2.5,3,4,5]/(1e6);


%DEBUG
% iterParamsCharges.Ek=[10e3, 100e3];
% iterParamsDevice.globalVa = [15,100];
% params.M                         = 100;
% params.N                         = 25;
% params.rPts                      = 200;
% params.zPts                      = 200;
% params.numOfParticles            = 10;


pcNames = fieldnames(iterParamsCharges);
pdNames = fieldnames(iterParamsDevice);


if isdir(['./simulations/',params.simGlobalName])
    rmdir(['./simulations/',params.simGlobalName],'s');
end
 
if (~exist(['./simulations/',params.simGlobalName,'/simulationsSummary'], 'dir'))
    mkdir (['./simulations/',params.simGlobalName,'/simulationsSummary'])
end

if(exist('SimulationsLog.txt', 'file'))
    delete('SimulationsLog.txt')
end

saveDefaultVals(params);
log = fopen('SimulationsLog.txt', 'wt');
fprintf(log, "Beginning %s... Time: %s \n",params.simGlobalName, datetime('now'));


results = struct();
results.in.globalDefaultParams = params;

%Writing the input values in the results vector
for pdNameInd = 1:numel(pdNames)
    name = pdNames{pdNameInd};
    results.in.device.(name) = iterParamsDevice.(name);
end

for pcNameInd = 1:numel(pcNames)
    name = pcNames{pcNameInd};
    results.in.charges.(name) = iterParamsCharges.(name);
end

for pdNameInd = 1:numel(pdNames)
    pdName = pdNames{pdNameInd};
    pdVec = iterParamsDevice.(pdName);
    for pdValInd = 1:length(pdVec)
        pdVal = pdVec(pdValInd);
        for pcNameInd = 1:numel(pcNames)
            pcName = pcNames{pcNameInd};
            pcVec = iterParamsCharges.(pcName);
            for pcValInd = 1:length(pcVec)  
                pcVal = pcVec(pcValInd);
                runParams = params;
                runParams.(pdName) = pdVal;
                runParams.(pcName) = pcVal;
                runParams.simName = sprintf('%s-[%d]-%s-[%d]',  ...
                                     pdName, pdVal, pcName, pcVal);
                runParams.simPdName = sprintf('%s-[%d]', pdName, pdVal);
                
                [focused_particles_percent, entryEmittance, exitEmittance, ...
                    randomSeed] = runSim(runParams);
                results.out.(pdName).(pcName).focused(pdValInd,pcValInd) = focused_particles_percent;
                results.out.(pdName).(pcName).entryEmittance(pdValInd,pcValInd) = entryEmittance;
                results.out.(pdName).(pcName).exitEmittance(pdValInd,pcValInd) = exitEmittance;
                results.out.(pdName).(pcName).randomSeed(pdValInd,pcValInd) = randomSeed;
                fprintf(log, "%s DONE, Time: %s \n",runParams.simName, datetime('now'));
            end
        end
    end
end

    
    for pdNameInd = 1:numel(pdNames)
        pdName = pdNames{pdNameInd};
        for pcNameInd = 1:numel(pcNames)
            pcName = pcNames{pcNameInd};
            fig = figure;
            entryEmittance = results.out.(pdName).(pcName).entryEmittance;
            exitEmittance = results.out.(pdName).(pcName).exitEmittance;
            deltaEmittance = 100*(exitEmittance-entryEmittance)./entryEmittance;
            plot(results.in.device.(pdName), deltaEmittance, '-o');
            tit = ['\DeltaEmittance vs. ', pdName, ' and ',  pcName];
            title(tit);
            xlabel(pdName);
            ylabel('\DeltaEmittance [%]');
            hold on
            
           %Plots the focused particles. Needs to be restored when we get a
           %closed formula for focused particles
           %TODO: plot it on two different graphs
%             plot(results.in.charges.(pcName), results.out.(pdName).(pcName).focused, '-o');
%             plot(results.in.device.(pdName), results.out.(pdName).(pcName).focused, '-o');
%             tit = ['Focused Particles vs. ', pdName, ' and ',  pcName];
%             title(tit);
%             xlabel(pdName);
%             ylabel('Focused particles [%]');
%             ylim([0,100]);
            
            
            legstr = string(results.in.charges.(pcName));
            legstr = strcat(strjoin([pcName, " = "]), legstr); 
            legend(legstr);
            savefig(fig, ['./simulations/', params.simGlobalName,'/simulationsSummary/', tit, '.fig']);
            close(fig);
        end
    end
    
    
    save(['./simulations/', params.simGlobalName,'/simulationsSummary/results.mat'], 'results');
    
    fprintf(log, "Simulations FINISHED, Time: %s \n", datetime('now'));
    fclose(log);
