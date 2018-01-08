function [ ] = fullSimEmittancePlot(simName, pdNames, pcNames, pdCombMat, pcCombMat, pdValMat, pcValMat)
pdCombOffset = 0;
globalSimPath = [ './simulations/',simName];
k=1;
for pd = 1:numel(pdNames)
    curPdNumOfVals = length(pdValMat.(pdNames{pd}));
    for pcCombIdx = 1:size(pcCombMat,2)
        curPcVals = pcCombMat(:,pcCombIdx);
        for combInnerOffset = 1:curPdNumOfVals
            curPdCombValIdx = pdCombOffset + combInnerOffset;
            %------ Load Device Results ------%
            curPdVals = pdCombMat(:,curPdCombValIdx);
            deviceResultsPath = [globalSimPath, '/DeviceResults/'];
            deviceValsStr = [];
            for i = 1:length(curPdVals)
                deviceValsStr = [deviceValsStr, pdNames{i}, ' = ',num2str(curPdVals(i)),', '];
            end
            deviceResultsPath = [deviceResultsPath, deviceValsStr(1:end-2), '.mat'];
%             load(deviceResultsPath);
            %----Load Trajectory Results -----%
            trajResPath = [globalSimPath, '/', deviceValsStr(1:end-2), ' - '];
            pcValsPath = [];
            for i = 1:length(curPcVals)
                pcValsPath = [pcValsPath, pcNames{i}, ' = ',num2str(curPcVals(i)),', '];
            end
            finalStrings{k} = [trajResPath, pcValsPath(1:end-2),'/ParticleTrajectory.mat'];
%             load(trajResPath);
            %------Plotting Current Emittance vs. Z-----%
            k =k+1;
        end
    end
    pdCombOffset = pdCombOffset+curPdNumOfVals;
end
            
end

