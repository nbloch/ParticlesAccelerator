function [ fig ] = plotEmittance( Z, X, Vz, Vx, Vy, zGrid)

rows = 1:size(X,1); 
i =1;
zAxis = 1:5:length(zGrid(1,:));
emittance = zeros(size(zAxis));
err = zeros(size(zAxis));
for j=zAxis
    [~, cols]= min(abs(Z-zGrid(1,j)),[],2);
    idxs = sub2ind(size(X), rows, cols');
    [BetaR, BetaZ, BetaY, ~, Gamma ] = getBG( Vz, Vx, Vy, cols');
    px = Gamma.*BetaR;
    pz = Gamma.*BetaZ;
    py = Gamma.*BetaY;
    [emittance(i), err(i)] = getEmittance( X(idxs), px, pz, py, Gamma);
    i=i+1;
end

fig = figure();
subplot(1,2,1)
plot(zGrid(1,1:5:end)*1e6, emittance);
title('Beams Emittance vs. Z Axis Position')
xlabel ('Position in Lens[\mum]')
ylabel ('Emittance[m x rad]')
subplot(1,2,2)
%For logarithmic scale: change "plot" below to "semilogy" and remove the
%multipliclation by 100 of Err in getEmittance
plot(zGrid(1,1:5:end)*1e6, err);
title('Beams Numeric Error')
xlabel ('Position in Lens[\mum]')
ylabel ('Numeric Error[%]')

end

