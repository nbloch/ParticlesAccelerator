function [focused_particles_percent, randomSeed] = runSim(params)

    %Loading all fields from params into the workspace
    save('params-tmp.mat', '-struct', 'params');
    load('params-tmp.mat');
    delete('params-tmp.mat');
    clear params;
    % END
    addpath(genpath([pwd]));
    
    if(exist('globalVa', 'var'))
        VaLeft = -globalVa;
        VaRight = globalVa;
    end
    
    if(exist('globalElectrodeRadius', 'var'))
        leftElectrodeRadius = globalElectrodeRadius;
        rightElectrodeRadius = globalElectrodeRadius;
    end

    deviceParams = {dimensionsOptimization, ...
              VaLeft, VaRight, electrodeWidth,leftElectrodeRadius,rightElectrodeRadius, deviceRadius, ...
              distanceBetweenElectrodes, use_bessel, repetitions, sideOffset,  ...
              M, N, rPts, zPts, convergeTh, growthTh};

%     filename = ['simulations/', simName, '/', 'deviceParams.mat'];
    filename = sprintf('./simulations/%s/deviceSimResults.mat', simGlobalName);
    calcPotential = ((exist (string(filename), 'file') == 0) && eraseOldSim);
    if(~calcPotential)
        load(filename, 'oldDeviceParams', 'V', 'zGrid', 'rGrid', 'Rq', 'Zq', 'Rb', 'Zb', 'MSE', 'Mbleft', 'Mbright');

         if(~isequal(oldDeviceParams, deviceParams))
            calcPotential = true;
        end
    end

    if(calcPotential)
        V =  [];
        zGrid =  [];
        rGrid =  []; 
        Rq = []; 
        Zq = [];
        Rb = []; 
        Zb = [];
        MSE=[];
        Mbleft = [];
        Mbright = [];
        disp('The potential is being computed...');
    else
        disp('Using precomputed potential');
    end

    if(~genNewSeed)
        if(~exist('randomSeed.mat', 'file'))
            error('The random seed file does not exist. Check the "Generate new randow seed checkbox.');
        end
        load('randomSeed.mat', 'randomSeed');
        rng(randomSeed);
    else
        randomSeed=rng;
        save('randomSeed.mat', 'randomSeed');
    end
    
    
    parsStr=   {'dimensionsOptimization', 'simulateTrajectory', 'BunchOrSingle', 'useAngle', 'simulateElectron', ...
              'VaLeft', 'VaRight', 'electrodeWidth','leftElectrodeRadius','rightElectrodeRadius', 'deviceRadius', ...
              'distanceBetweenElectrodes', 'use_bessel', 'repetitions', 'sideOffset',  ...
              'M', 'N', 'rPts', 'zPts', 'convergeTh', 'growthTh', ...
              'entryVel', 'entryAngle', 'entryR', 'q', 'm', 'axialEntryVelocity', 'radialEntryVelocity', 'useRelativeTrajectory', ...
              'electrodeProximityThresh', 'numOfParticles', 'lowAxialVel', 'highAxialVel',...
              'lowRadialVel', 'highRadialVel','simulatePhaseSpace', ...
              'exitRthresh', 'calcPotential', 'V', 'zGrid', 'rGrid', ...
              'Rq', 'Zq', 'Rb', 'Zb', 'MSE', 'Mbleft', 'Mbright','recordPhaseSpace',...
              'beamInitialRadius'};
         
    for i=1:length(parsStr(:))
        varname = parsStr{i};
        simPars.(varname)=eval(varname);
    end
          
          
          
          
    [ V, X, Y, U, W, MSE, fig, zGrid, rGrid, Rq, Zq, Rb, Zb, Ez, Er, Mbleft, Mbright, focused_particles ] = FullEinzelSim(simPars);

     simDir = ['simulations/', simGlobalName,'/', simName];
     if ~isdir(simDir)
         mkdir(simDir);
     end
     
     cd(simDir);
     
     if(~isdir('../DeviceResults'))
         mkdir('../DeviceResults')
     end
     

     if (calcPotential) 
         oldDeviceParams = deviceParams;
         save('../deviceSimResults.mat', 'oldDeviceParams',  'V','MSE',...
             'zGrid', 'rGrid', 'Rq', 'Zq', 'Rb', 'Zb', 'Mbleft',...
             'Mbright', 'Ez', 'Er');
         if (SingleSim)
             if(dimensionsOptimization)
                 savefig(fig(1),'MSEFigure.fig','compact');
             end
             savefig(fig(2),'ChargeDistribution.fig','compact');
         else
             save(['../DeviceResults/', simPdName,'.mat'],  ...
             'V','MSE', 'zGrid', 'rGrid', 'Rq', 'Zq', 'Rb', 'Zb',...
             'Mbleft', 'Mbright', 'Ez', 'Er');
         end
     end 
     
     if(SingleSim)
         if(simulatePhaseSpace) 
             savefig(fig(3),'phaseSpace.fig','compact');
         end
         savefig(fig(4),'FullSim.fig','compact');
     end
     
     save('ParticleTrajectory.mat', 'X', 'Y', 'U', 'W')
    

%      if (strcmp(class(probFig), 'matlab.graphics.Graphics')) %in case there are no problematic particles.
     if isdir('./problematicParticles/')
        rmdir('./problematicParticles','s');
     end
     if isdir('../../../problematicParticles/')
         movefile('../../../problematicParticles/')
     end  
     
     if (recordPhaseSpace)
        movefile('../../../phaseSpace.avi')
     end
     
     cd ../../..

    focused_particles_percent = 100*focused_particles/numOfParticles;

end

