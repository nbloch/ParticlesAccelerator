function [ N, M, MSE, fig ] = FindOperatingPoint(VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                       deviceRadius, distanceBetweenElectrodes,repetitions, convergeTh, growthTh )
                                       
%Gets the geometry of the problem and the number of boundaries required,
%and finds the minimal number of charges that minimizes the MSE
Nstep = 10;
Mstep = 10;

minNnum = 200;
maxNnum = 1000;

maxMperN = 1000;
N_vec = minNnum : Nstep : maxNnum;

MSE   = NaN*ones(length(N_vec), maxMperN);
M_mat = NaN*ones(length(N_vec), maxMperN);

%------------------------------------
%Calculate MSE for each configuration
%------------------------------------
for i=1:length(N_vec)
    tmp   = NaN*ones(1, maxMperN);
    M_vec = NaN*ones(1, maxMperN);
    M = 2*N_vec(i);
    for j = 1:maxMperN
        [~, ~, ~, ~, ~, tmp(j)] = getChargesRepetitive(VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                              deviceRadius, distanceBetweenElectrodes, N_vec(i), M, repetitions);                             

        M_vec(j) = M;
        M = M + Mstep;
        if ((j >= 2) && ( ((((tmp(j-1) - tmp(j))> 0) && (tmp(j-1) - tmp(j))< convergeTh)) || ((tmp(j) - tmp(j-1)) > growthTh))) 
            break
        end
    end
    M_mat(i,:) = M_vec;
    MSE(i,:) = tmp;
end


%------------------------------
%Plotting
%------------------------------
for i = 1:length(N_vec)
    legstr{i} = sprintf('N = %s',num2str(N_vec(i)));
end

save('MSEDATA.mat', 'MSE','M_mat', 'N_vec');
fig = figure();
    for i = 1:length(N_vec)
        plot(M_mat(i,:), MSE(i,:), '-o');
        
        hold on;
    end
    hold off;
    title(sprintf('Dimensions Optimization Algorithm For %d Unit Cells', repetitions));
    xlabel('M- Num of Boundaries Conditions');
    ylabel('MSE [%]');
    hleg = legend(legstr,'Location', 'best');
    htitle = get(hleg, 'Title');
    set(htitle,'String', 'Number Of Charges Per Cell');


%------------------------------
%Choose Minimum Value
%------------------------------    
min_val = min(min(MSE));
[Nind, Mind] = find(MSE == min_val);
N = N_vec(Nind);
M = M_mat(Nind,Mind);

end

