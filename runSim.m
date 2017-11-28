function [focused_particles_percent, randomSeed] = runSim(params)

    %Loading all fields from params into the workspace
    save('params-tmp.mat', '-struct', 'params');
    load('params-tmp.mat');
    delete('params-tmp.mat');
    clear params;
    % END
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
    filename = sprintf('./simulations/%s/deviceParams.mat', simGlobalName);
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
    end
    if(calcPotential)
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
              'Rq', 'Zq', 'Rb', 'Zb', 'MSE', 'Mbleft', 'Mbright','recordPhaseSpace'};
         
    for i=1:length(parsStr(:))
        varname = parsStr{i};
        simPars.(varname)=eval(varname);
    end
          
          
          
          
    [ V, X, Y, U, W, MSE, fig, zGrid, rGrid, Rq, Zq, Rb, Zb, Ez, Er, probFig, Mbleft, Mbright, focused_particles ] = FullEinzelSim(simPars);

     simDir = ['simulations/', simGlobalName,'/', simName];
     if ~isdir(simDir)
         mkdir(simDir);
     end
     
     cd(simDir);


     if (calcPotential)
%          folderExist = isdir (['./','simulations']);
%          if (~folderExist)
%              mkdir simulations
%          end
%          cd 'simulations'
%          if(~isdir(['./',simGlobalName]))
%              mkdir(simGlobalName);
%          end
% 
%          cd(simGlobalName);
         
         oldDeviceParams = deviceParams;
         save('../deviceParams.mat', 'oldDeviceParams',  'V', 'X', 'Y','U', 'W','MSE', 'zGrid', 'rGrid', 'Rq', 'Zq', 'Rb', 'Zb', 'Mbleft', 'Mbright', 'Ez', 'Er');

%          if (isdir(['./',simName]))
%              rmdir(['./',simName], 's');
%          end
%          
% 
%          mkdir(simName);
%          cd(simName);
         

         if(dimensionsOptimization)
             savefig(fig(1),'MSEFigure.fig','compact');
         end
         savefig(fig(2),'ChargeDistribution.fig','compact');
     end 

     if(simulatePhaseSpace) 
         savefig(fig(3),'phaseSpace.fig','compact');
     end
     savefig(fig(4),'FullSim.fig','compact');

%      if (strcmp(class(probFig), 'matlab.graphics.Graphics')) %in case there are no problematic particles.
    if (strcmp(class(probFig), 'matlab.ui.Figure')) %in case there are no problematic particles.
         for i = 1:length(probFig)
             savefig(probFig(i),sprintf('ProblematicParticle-figNum-%d', i),'compact');
         end
     end
     
     if (recordPhaseSpace)
        movefile('../../../phaseSpace.avi')
     end
     
     cd ../../..

    focused_particles_percent = 100*focused_particles/numOfParticles;

end

