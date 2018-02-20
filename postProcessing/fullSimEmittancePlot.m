function [ ] = fullSimEmittancePlot(simName, pdNames, pcNames, pdCombMat, pcCombMat, pdValMat, pcValMat)
pdCombOffset = 0;
globalSimPath = [ './simulations/',simName];
deviceResultsPath = [globalSimPath, '/DeviceResults/'];
k=1;
log = fopen('emittanceVsZlog.txt', 'wt');
fprintf(log, "Starting, Time: %s \n",datetime('now'));
emittanceLenFactor = 200000;
for pd = 1:numel(pdNames)
    curPdNumOfVals = length(pdValMat.(pdNames{pd}));
    for pcCombIdx = 1:size(pcCombMat,2)
        curPcVals = pcCombMat(:,pcCombIdx);
        pcValsPath = [];
        for i = 1:length(curPcVals)
            pcValsPath = [pcValsPath, pcNames{i}, ' = ',num2str(curPcVals(i)),', '];
        end
        electrodesZ = NaN*ones(19,curPdNumOfVals);
        electrodPosEmit = NaN*ones(19,curPdNumOfVals);
        emittance = NaN.*ones(emittanceLenFactor,curPdNumOfVals);
        zPtsCalc = NaN.*ones(emittanceLenFactor,curPdNumOfVals);
        legStr = [];
        for combInnerOffset = 1:curPdNumOfVals
            curPdCombValIdx = pdCombOffset + combInnerOffset;
            %------ Load Device Results ------%
            curPdVals = pdCombMat(:,curPdCombValIdx);
            deviceValsStr = [];
            for i = 1:length(curPdVals)
                deviceValsStr = [deviceValsStr, pdNames{i}, ' = ',num2str(curPdVals(i)),', '];
            end
            deviceResultsMatPath = [deviceResultsPath, deviceValsStr(1:end-2), '.mat'];
            load(deviceResultsMatPath);
            
            %----Load Trajectory Results -----%
            trajResPath = [globalSimPath, '/', deviceValsStr(1:end-2), ' - '];
            finalStrings{k} = [trajResPath, pcValsPath(1:end-2),'/ParticleTrajectory.mat'];
            fprintf(log, "%s \n", finalStrings{k});
            load( finalStrings{k});
            
            %------Calculating Current Emittance vs. Z-----%
            [ curEmittance, curzPtsCalc, ~] = calcEmiitanceVsZ( Z, X, Vz, Vx, Vy, zGrid );
            curLen = length(curEmittance);
            emittance(1:curLen,combInnerOffset) = curEmittance;
            zPtsCalc(1:curLen, combInnerOffset) = curzPtsCalc;
            
            %------Calculating Points For Electrodes Position Plot----%
            curElecZ = unique(Zq(1,:));
            numOfEle = length(curElecZ);
            electrodesZ(1:numOfEle,combInnerOffset) = curElecZ;
             
            [~, lensZidx] = min(abs(repmat(curElecZ',1,curLen) - repmat(curzPtsCalc,numOfEle,1)),[],2);
            electrodPosEmit(1:numOfEle,combInnerOffset) = curEmittance(lensZidx);
            legStr{combInnerOffset} = sprintf('%s = %s',pdNames{pd},num2str(curPdVals(pd)) );
            k =k+1;
        end
         %----------Centering around First Electrode --------%
        leftestElec = min(electrodesZ(1,:));
        perElecShiftFac = abs(electrodesZ(1,:) - leftestElec);
        
        elecZAlignMat = repmat(perElecShiftFac,size(electrodesZ,1),1);
        zPtsAlignMat  = repmat(perElecShiftFac,size(zPtsCalc,1),1);   
  
        electrodesZ = electrodesZ - elecZAlignMat;
        zPtsCalc = zPtsCalc -zPtsAlignMat;
        
        positiveElectrodePos  = electrodesZ(1:2:end,:)*1e6;
        positiveElectrodeEmit = electrodPosEmit(1:2:end,:);
        negativeElectrodePos  = electrodesZ(2:2:end,:)*1e6;
        negativeElectrodeEmit = electrodPosEmit(2:2:end,:);
        legStr{combInnerOffset+1} = sprintf('Positive Electrode Position');
        legStr{combInnerOffset+2} = sprintf('Negative Electrode Position');
        
         %----------Plotting --------%
        fig = figure();
        plot(zPtsCalc*1e6, emittance);
        hold on
        plot(positiveElectrodePos(:), positiveElectrodeEmit(:),'r+')
        plot(negativeElectrodePos(:), negativeElectrodeEmit(:),'bo')
        title(['Emittance Vs. Z - Changing Lens Parameter: ', pdNames(pd)]);
        xlabel('Z [\mum]')
        ylabel('\epsilon(z)[m x rad]')
        legend(legStr,'Location', 'northeast')
        
        %----------Saving Figure --------%
        cd ([globalSimPath,'/EmittanceVsZ-Summary'])
        otherPdVals = [];
        for i = 1:length(curPdVals)
            if(i == pd)
                continue;
            else
            otherPdVals = [otherPdVals, pdNames{i}, ' = ',num2str(curPdVals(i)),', '];
            end
        end
        figName = ['Tested Param- ', pdNames{pd},', Othe Parms Vals- ' ,otherPdVals, pcValsPath,'EmittanceVsZ.fig'];
        savefig(fig, figName,'compact');
        cd ../../..
        close(fig)
    end
    pdCombOffset = pdCombOffset+curPdNumOfVals;
end
fprintf(log, "Simulations FINISHED, Time: %s \n", datetime('now'));
fclose(log);   
end

