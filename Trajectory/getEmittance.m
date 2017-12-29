function  emittance = getEmittance( x, px, gamma )
xMean = mean(x);
pxMean = mean(px);
emittance = 4*sqrt(mean((x-xMean).^2)*mean(((px-pxMean)./gamma).^2) ...
                            - mean(((x-xMean).*(px-pxMean)./gamma)).^2);
end

