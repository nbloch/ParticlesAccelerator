function [pMat, pMatNames] = buildParamMat(iterParams, params, combination)
    %pMatNames holds the name of the parameter we iterate on (only relevant
    %for combination = true

    pNames = fieldnames(iterParams);
    pMatNames = [];
    %pcList is a cell array that stores a list of vectors that contains the
    %possible values for each iteration parameter
    pList = cell(1,numel(pNames));
    totalParamValues = 0;
    for pNameInd = 1:numel(pNames)
        pName = pNames{pNameInd};
        pVec = iterParams.(pName);
        pList{pNameInd} = pVec;
        totalParamValues = totalParamValues + length(pVec);
    end

    %pcMat initialization based on pcList
    if combination
        pMat = combvec(pList{:});
    else
        i = 1;
        pMat = zeros(numel(pList), totalParamValues);
        pMatNames = strings(1, totalParamValues);
        for pNameInd = 1:numel(pNames)
            pName = string(pNames{pNameInd});
            pVec = iterParams.(pName);
            % We first initialize the whole line whith the default param
            pMat(pNameInd, :) = params.(pName);
            % Then we replace the relevant elements in the matrix
            pMat(pNameInd, i:i+numel(pVec)-1) = pVec;
            pMatNames(i:i+numel(pVec)-1) = pName;
            i = i + numel(pVec);
        end
    end

end

