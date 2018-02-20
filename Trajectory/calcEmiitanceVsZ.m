function [ emittance, zPtsCalc, err] = calcEmiitanceVsZ( Z, X, Vz, Vx, Vy, zGrid )
%CALCEMIITANCEVSZ Summary of this function goes here
%   Detailed explanation goes here
rows = 1:size(X,1); 
% i =1;
zAxisIdx = 1:length(zGrid(1,:));
emittance = zeros(size(zAxisIdx));
err = zeros(size(zAxisIdx));
zPtsCalc = zGrid(1,zAxisIdx);
parfor j=zAxisIdx
    [~, cols]= min(abs(Z-zGrid(1,j)),[],2);
    idxs = sub2ind(size(X), rows, cols');
    [BetaR, BetaZ, BetaY, ~, Gamma ] = getBG( Vz, Vx, Vy, cols');
    px = Gamma.*BetaR;
    pz = Gamma.*BetaZ;
    py = Gamma.*BetaY;
    [tmpEmit, tmpErr] = getEmittance( X(idxs), px, pz, py, Gamma);
    err(j) = tmpErr;
    emittance(j) = tmpEmit;
end


end

