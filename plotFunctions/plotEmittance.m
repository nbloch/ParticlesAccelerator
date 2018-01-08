function [ fig ] = plotEmittance( zPtsCalc, emittance, err )

fig = figure();
subplot(1,2,1)
plot(zPtsCalc*1e6, emittance);
title('Beams Emittance vs. Z Axis Position')
xlabel ('Position in Lens[\mum]')
ylabel ('Emittance[m x rad]')
subplot(1,2,2)
%For logarithmic scale: change "plot" below to "semilogy" and remove the
%multipliclation by 100 of Err in getEmittance
plot(zPtsCalc*1e6, err);
title('Beams Numeric Error')
xlabel ('Position in Lens[\mum]')
ylabel ('Numeric Error[%]')

end

