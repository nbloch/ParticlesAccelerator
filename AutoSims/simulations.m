%"results" struct documentation:
% 2 fields:
% -in:
%   3 fields:
%   -globalDefaultParams: The default parameters of the simulation (not the
%   iteration parameters)
%   - device: device iteration parameters
%   - charges: charges iteration parameters
%   Both of the above fields contain the names of the parameters we iterate
%   on along with their values (stored as vectors)
% -out:
%   N fields (N = num of different device parameters we simulate. e.g.
%   param1 = [1 2], param2 = [5 6 7] --> N = 2). Each field is a device
%   parameter name. In each of them:
%       M fields (M = how many SETS of charge parameter are simulated. e.g.
%       param1 = [1 2 3], param2 = [5 6 7], combinatory --> N = 9, not
%       combinatory --> N = 6. Each fields name contains the whole set of
%       charge params, separated with X and underscores (special characters
%       are not allowed in field names...). 
%       Example: 
%       If we simulate param1= [1.5 2], param2 = [3e5], we'll have two
%       fields whose names will be: param1X1_5Xparam2X3e5,
%       param1X2Xparam2X3e5
%       Each of the above fields contains the results of each simulation


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
params.globalVa                  = 15;

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
params.plotFullSim              = false;
params.plotPhaseSpace           = false;
params.plotPhaseSpaceVideo      = false;
params.plotChargeDistribution   = false;
params.plotProblematicParticles = false;
params.plotEmittanceVsZ         = false;

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
params.genNewSeed    = false;
params.savePlots     = true;
%params.simPdName     = 'Repetitions';
params.SingleSim     = false;

%How do we combine the pc and pd parameters in the simulations? 
%True = All cominations 
%False = each parameter solely with other parameters default values

params.pcCombinations=true;
%Plotting of pdCombinations not implemeted yet
params.pdCombinations=false;



if params.pdCombinations == true
    error('Combinatorial device parameters not implemented yet');
end

paramsFields = fieldnames(params);

iterParamsCharges = struct();
iterParamsDevice = struct();

% %charge parameters iteration variables definitions (must be declared as ROW
% %vectors)
% 
% iterParamsCharges.Ek                        = [10e3, 100e3];    %eV units
% iterParamsCharges.maxInitMoment             = [1e-3, 1e-2];
% 
% 
% %device parameters iteration variables definitions(must be declared as
% %ROW vectors)
% 
% iterParamsDevice.globalVa = [10,15,30,50,100,300,500,750,1000,2000];
% %iterParamsDevice.repetitions = [1,2,3,5,7,9];
% % iterParamsDevice.distanceBetweenElectrodes=[0.5,0.75,1,1.25,2]/1e6;
% % iterParamsDevice.globalElectrodeRadius=[0.25,0.5,0.75,1,1.5,2]/1e6;
% % iterParamsDevice.lensPreOffset = [1,1.5,2,2.5,3,4,5]/(1e6);


% DEBUG
iterParamsCharges.Ek=[10e3, 100e3, 200e3];
iterParamsCharges.maxInitMoment     = [1e-3, 1e-2];
iterParamsDevice.globalVa = [15,100];
iterParamsDevice.repetitions = [1,2];
%iterParamsDevice.lensPreOffset = [1,1.5]/(1e6);
params.M                         = 100;
params.N                         = 25;
params.rPts                      = 200;
params.zPts                      = 200;
params.numOfParticles            = 10;


% Cell arrays that contain the names of the iteration parameters (both for
% charges and device)
pcNames = fieldnames(iterParamsCharges);
pdNames = fieldnames(iterParamsDevice);
params.pdNames = pdNames;
params.pcNames = pcNames;

%Contains the values for the iteration parameters. Each column correspond to a new
%simulation, each row to a parameter name
[pcMat, pcMatNames] = buildParamMat(iterParamsCharges, params, params.pcCombinations);
[pdMat, pdMatNames] = buildParamMat(iterParamsDevice, params, params.pdCombinations);

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

% %OLD IMPLEMENTATION (DOESN'T SUPPORT COMBINATORY PARAMETERS)
% for pdNameInd = 1:numel(pdNames)
%     pdName = pdNames{pdNameInd};
%     pdVec = iterParamsDevice.(pdName);
%     for pdValInd = 1:length(pdVec)
%         pdVal = pdVec(pdValInd);
%         for pcNameInd = 1:numel(pcNames)
%             pcName = pcNames{pcNameInd};
%             pcVec = iterParamsCharges.(pcName);
%             for pcValInd = 1:length(pcVec)  
%                 pcVal = pcVec(pcValInd);
%                 runParams = params;
%                 runParams.(pdName) = pdVal;
%                 runParams.(pcName) = pcVal;
%                 runParams.simName = sprintf('%s-[%d]-%s-[%d]',  ...
%                                      pdName, pdVal, pcName, pcVal);
%                 runParams.simPdName = sprintf('%s-[%d]', pdName, pdVal);
%                 
%                 [focused_particles_percent, entryEmittance, exitEmittance, ...
%                     randomSeed] = runSim(runParams);
%                 results.out.(pdName).(pcName).focused(pdValInd,pcValInd) = focused_particles_percent;
%                 results.out.(pdName).(pcName).entryEmittance(pdValInd,pcValInd) = entryEmittance;
%                 results.out.(pdName).(pcName).exitEmittance(pdValInd,pcValInd) = exitEmittance;
%                 results.out.(pdName).(pcName).randomSeed(pdValInd,pcValInd) = randomSeed;
%                 fprintf(log, "%s DONE, Time: %s \n",runParams.simName, datetime('now'));
%             end
%         end
%     end
% end


runParams = params;
%Each column of pcMat and of pdMat represent a set of pc and pd parameters
%respectively and thus represent one simulation

%Iteration over all pdMat columns
for pdSetInd = 1 : size(pdMat, 2)
    pdSet = pdMat(:, pdSetInd);
    %Setting the new input parameters in runParams
    for pdNameInd = 1:numel(pdSet)
        pdName = pdNames{pdNameInd};
        runParams.(pdName) = pdSet(pdNameInd);
    end
    %Iteration over all pcMat columns
    for pcSetInd = 1 : size(pcMat, 2)
        pcSet = pcMat(:, pcSetInd);
        %Setting the new input parameters in runParams
        for pcNameInd = 1:numel(pcSet)
            pcName = pcNames{pcNameInd};
            runParams.(pcName) = pcSet(pcNameInd);
        end

        pdTitle = buildParamString(pdNames, pdSet, false);
        pcTitle = buildParamString(pcNames, pcSet, false);
    
        runParams.simPdName = char(pdTitle);        
        runParams.simPcName = char(pcTitle);   %Not in use for now
        runParams.simName = char(pdTitle + " - " + pcTitle);
        
        %Will define the field name in results
        pdName = char(pdMatNames(pdSetInd));
        pcString = char(buildParamString(pcNames, pcSet, true));
        
        
        [focused_particles_percent, entryEmittance, exitEmittance, ...
            randomSeed] = runSim(runParams);
        
        %The condition is for 'focused' but can be generalized to the other
        %parameters
        if ~isfield(results, 'out') || ...
                ~isfield(results.out, pdName) ||...
                ~isfield(results.out.(pdName), pcString) || ...
                ~isfield(results.out.(pdName).(pcString), 'focused')
            
            results.out.(pdName).(pcString).focused = [];
            results.out.(pdName).(pcString).entryEmittance = [];
            results.out.(pdName).(pcString).exitEmittance = [];
            results.out.(pdName).(pcString).randomSeed = [];
        end
        
        results.out.(pdName).(pcString).focused = ...
            [results.out.(pdName).(pcString).focused, focused_particles_percent];
        results.out.(pdName).(pcString).entryEmittance = ...
            [results.out.(pdName).(pcString).entryEmittance, entryEmittance];
        results.out.(pdName).(pcString).exitEmittance = ...
            [results.out.(pdName).(pcString).exitEmittance, exitEmittance];
        results.out.(pdName).(pcString).randomSeed = ...
            [results.out.(pdName).(pcString).randomSeed, randomSeed];
        results.out.(pdName).(pcString).pcTitle = pcTitle;
        
        fprintf(log, "%s DONE, Time: %s \n",runParams.simName, datetime('now'));
    end
end

    for pdNameInd = 1:numel(pdNames)
        pdName = pdNames{pdNameInd};
        pcStrings = fieldnames(results.out.(pdName));
        
        fig = figure;
        tit = ['\DeltaEmittance vs. ', pdName];
        title(tit);
        xlabel(pdName);
        ylabel('\DeltaEmittance [%]');
        hold on
        legstr = cell(1, numel(pcStrings));
        
        for pcNameInd = 1:numel(pcStrings)
            pcString = pcStrings{pcNameInd};
            entryEmittance = results.out.(pdName).(pcString).entryEmittance;
            exitEmittance = results.out.(pdName).(pcString).exitEmittance;
            deltaEmittance = 100*(exitEmittance-entryEmittance)./entryEmittance;
            plot(results.in.device.(pdName), deltaEmittance, '-o');

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
             
            legstr{pcNameInd} = char(results.out.(pdName).(pcString).pcTitle);
            
        end
        
        legend(legstr);
        savefig(fig, ['./simulations/', params.simGlobalName,'/simulationsSummary/', tit, '.fig']);
        close(fig);
    end
    
    
    save(['./simulations/', params.simGlobalName,'/simulationsSummary/results.mat'], 'results');
    
    fprintf(log, "Simulations FINISHED, Time: %s \n", datetime('now'));
    fclose(log);
    
    [ ~ ] = fullSimEmittancePlot(params.simGlobalName, pdNames, pcNames, pdMat, pcMat, iterParamsDevice, iterParamsCharges);
