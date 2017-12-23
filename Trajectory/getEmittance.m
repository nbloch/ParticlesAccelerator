function  emittance = getEmittance( x, px, gamma, lambda0 )
xMean = mean(x);
pxMean = mean(px);
emittance = 4*lambda0*sqrt(mean((x-xMean).^2)*mean(((px-pxMean)./gamma).^2) ...
                            - mean(((x-xMean).*(px-pxMean)./gamma)).^2);
end

