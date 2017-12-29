function [ fig ] = plotEmittance( Z, X, Vz, Vx, zGrid)

rows = 1:size(X,1); 
i =1;
zAxis = 1:5:length(zGrid(1,:));
emittance = zeros(size(zAxis));
for j=zAxis
    [~, cols]= min(abs(Z-zGrid(1,j)),[],2);
    idxs = sub2ind(size(X), rows, cols');
    [BetaR, ~, ~, Gamma ] = getBG( Vz, Vx, cols');
    px = Gamma.*BetaR;
    emittance(i) = getEmittance( X(idxs), px, Gamma );
    i=i+1;
end

fig = figure();
plot(zGrid(1,1:5:end)*1e6, emittance)
title('Beams Emittance vs. Z Axis Position')
xlabel ('Position in Lens[\mum]')
ylabel ('Emittance')

end
