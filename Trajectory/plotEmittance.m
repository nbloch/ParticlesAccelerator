function [ fig ] = plotEmittance( Z, X, Vz, Vx, zGrid)

rows = 1:size(X,1); 
i =1;
zAxis = 1:5:length(zGrid(1,:));
emittance = zeros(size(zAxis));
for j=zAxis
    [~, cols]= min(abs(Z-zGrid(1,j)),[],2);
    idxs = sub2ind(size(X), rows, cols');
    [BetaR, BetaZ, ~, Gamma ] = getBG( Vz, Vx, cols');
    px = Gamma.*BetaR;
    pz = Gamma.*BetaZ;
    [emittance(i), err(i)] = getEmittance( X(idxs), px, pz, Gamma);
    i=i+1;
end

fig = figure();
subplot(1,2,1)
plot(zGrid(1,1:5:end)*1e6, emittance);
title('Beams Emittance vs. Z Axis Position')
xlabel ('Position in Lens[\mum]')
ylabel ('Emittance[\mum x radians]')
subplot(1,2,2)
plot(zGrid(1,1:5:end)*1e6, err);
title('Beams Numeric Error')
xlabel ('Position in Lens[\mum]')
ylabel ('Numeric Error[%]')

end

