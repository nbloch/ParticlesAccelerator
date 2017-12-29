function [focused_particles_percent, entryEmittance, exitEmittance, randomSeed] = runSim(params) %#ok<INUSD>
    %including all sub folders
    addpath(genpath(pwd));
    
    %-------------------------------------------------%
    %Loading all fields from params into the workspace
    %-------------------------------------------------%
    save('params-tmp.mat', '-struct', 'params');
    load('params-tmp.mat');
    delete('params-tmp.mat');
    clear params;
    
    %---------------------------------%
    %Auto Sims Parameters Modification:  
    %---------------------------------%
    if(exist('globalVa', 'var'))
        VaLeft = -globalVa;
        VaRight = globalVa;
    end
    
    if(exist('globalElectrodeRadius', 'var'))
        leftElectrodeRadius = globalElectrodeRadius;
        rightElectrodeRadius = globalElectrodeRadius;
    end
    %-------------------------------------%
    %Caching - Use Of Precomputed Potential
    %-------------------------------------%
    deviceParams = {dimensionsOptimization, ...
              VaLeft, VaRight, electrodeWidth,leftElectrodeRadius,rightElectrodeRadius, deviceRadius, ...
              distanceBetweenElectrodes, repetitions, lensPreOffset, lensPostOffset, M, N, rPts, zPts, convergeTh, growthTh};

    filename = sprintf('./simulations/%s/deviceSimResults.mat', simGlobalName);
    calcPotential = ((exist (string(filename), 'file') == 0) && eraseOldSim);
    if(~calcPotential)
        load(filename, 'oldDeviceParams', 'V', 'zGrid', 'rGrid', 'Rq', 'Zq', 'Rb', 'Zb', 'MSE', 'Mbleft', 'Mbright','NLeft', 'Q');
        if(~isequal(oldDeviceParams, deviceParams)) %#ok<NODEF>
            calcPotential = true;
        end
    end

    if(calcPotential)
        V =  []; %#ok<NASGU>
        zGrid =  []; %#ok<NASGU>
        rGrid =  [];  %#ok<NASGU>
        Rq = [];  %#ok<NASGU>
        Zq = []; %#ok<NASGU>
        Rb = [];  %#ok<NASGU>
        Zb = []; %#ok<NASGU>
        MSE=[]; %#ok<NASGU>
        Mbleft = []; %#ok<NASGU>
        Mbright = []; %#ok<NASGU>
        Q= []; %#ok<NASGU>
        NLeft = []; %#ok<NASGU>
        disp('The potential is being computed...');
    else
        disp('Using precomputed potential');
    end
    %---------------------------%
    %New Seed Randomization
    %---------------------------%
    if(~genNewSeed)
        if(~exist('randomSeed.mat', 'file'))
            error('The random seed file does not exist. Check the "Generate new randow seed checkbox.');
        end
        load('randomSeed.mat', 'randomSeed');
        rng(randomSeed); %#ok<NODEF>
    else
        randomSeed=rng;
        save('randomSeed.mat', 'randomSeed');
    end
    
    %---------------------------%
    %Call Simulator
    %---------------------------%
    parsStr ={... %device Parameters
              'deviceRadius', 'lensPreOffset', 'lensPostOffset', 'distanceBetweenElectrodes', 'leftElectrodeRadius','rightElectrodeRadius',...
              'electrodeWidth','repetitions', 'VaLeft', 'VaRight',... 
              ... %Resolution Parameters
              'M', 'N', 'rPts', 'zPts','dimensionsOptimization', 'convergeTh', 'growthTh', ...
              ... %what To Simulate
              'simulateSingleParticle', 'simulateMultiParticles', ...
              ... %Figures
              'plotFullSim', 'plotPhaseSpace', 'plotPhaseSpaceVideo', 'plotChargeDistribution', 'plotProblematicParticles',...
              'plotEmittanceVsZ',...
              ... %Particle Parameters
              'electrodeProximityThresh', 'q', 'm', 'simulateElectron', ...
              ... %Multi Particle
              'numOfParticles', 'beamInitialRadius', 'maxInitMoment', 'Ek',...
              ... %Single Particle
              'axialEntryVelocity', 'radialEntryVelocity', 'entryR','exitRthresh', ...
              ... %Cached Results
              'calcPotential', 'V', 'zGrid', 'rGrid', 'Rq', 'Zq', 'Rb', 'Zb', 'MSE', 'Mbleft', 'Mbright','NLeft','Q'};         
    for i=1:length(parsStr(:))
        varname = parsStr{i};
        simPars.(varname)=eval(varname);
    end

    [ V, Z, X, Vz, Vx, MSE, fig, zGrid, rGrid, Rq, Zq, Rb, Zb, Ez, Er, ...
        Mbleft, Mbright, NLeft, focused_particles, entryEmittance, ...
        exitEmittance, Q ] ...
        = FullEinzelSim(simPars); %#ok<ASGLU>

    %---------------------------%
    %Prepare File system
    %---------------------------%
     simDir = ['simulations/', simGlobalName,'/', simName];
     if ~isdir(simDir)
         mkdir(simDir);
     end
     
     cd(simDir);
     
     if(~isdir('../DeviceResults'))
         mkdir('../DeviceResults')
     end
     
     if isdir('./problematicParticles/')
        rmdir('./problematicParticles','s');
     end
     
    %---------------------------%
    %Save Figures and Data
    %---------------------------%
     if (calcPotential)
         oldDeviceParams = deviceParams; %#ok<NASGU>
         %for cahching
         save('../deviceSimResults.mat', 'oldDeviceParams',  'V','MSE',...
             'zGrid', 'rGrid', 'Rq', 'Zq', 'Rb', 'Zb', 'Mbleft',...
             'Mbright', 'Ez', 'Er', 'Q', 'NLeft');
         %saving lens results
         if (SingleSim) 
             movefile('../deviceSimResults.mat');
         else
             save(['../DeviceResults/', simPdName,'.mat'],  ...
                 'V','MSE', 'zGrid', 'rGrid', 'Rq', 'Zq', 'Rb', 'Zb',...
                 'Mbleft', 'Mbright', 'Ez', 'Er', 'Q', 'NLeft');
         end
     end 
     
     save('ParticleTrajectory.mat', 'Z', 'X', 'Vz', 'Vx')

     if (savePlots)
         if(plotFullSim);                                savefig(fig.FullSim, 'FullSim.fig','compact');                  end
         if(plotChargeDistribution);                     savefig(fig.ChargeDist, 'ChargeDistribution.fig','compact');    end
         if(plotPhaseSpace   & (fig.PhaseSpace  ~=0));   savefig(fig.PhaseSpace, 'PhaseSpace.fig','compact');            end
         if(plotEmittanceVsZ & (fig.EmittanceVsZ~=0));   savefig(fig.EmittanceVsZ, 'EmittanceVsZ.fig','compact');        end
         if(dimensionsOptimization);                     savefig(fig.MSECalc, 'MSECalc.fig','compact');                  end
         if(plotPhaseSpaceVideo & exist('../../../phaseSpace.avi','file'));        movefile('../../../phaseSpace.avi');        end
         if(plotProblematicParticles & isdir('../../../problematicParticles/'));   movefile('../../../problematicParticles/'); end 
     end

     cd ../../..
     
     focused_particles_percent = 100*focused_particles/numOfParticles;

end

