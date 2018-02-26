function [numErr, emitFig] = plotEmittance( zPtsCalc, emittance, err )
axisFont = 14;
emitFig = figure();
plot(zPtsCalc*1e6, emittance);
pbaspect([1 1 1]);
title('Beams Emittance vs. Z Axis Position')
ax = gca;
ax.TitleFontSizeMultiplier = 2;
xlabel ('Position in Lens[\mum]','FontSize',axisFont)
xlim([zPtsCalc(1)*1e6, zPtsCalc(end)*1e6]);
ylabel ('Emittance[m x rad]','FontSize',axisFont)
numErr = figure();
%For logarithmic scale: change "plot" below to "semilogy" and remove the
%multipliclation by 100 of Err in getEmittance
plot(zPtsCalc*1e6, err);
pbaspect([1 1 1]);
title('Beams Numeric Error')
ax = gca;
ax.TitleFontSizeMultiplier = 2;
xlabel ('Position in Lens[\mum]','FontSize',axisFont)
xlim([zPtsCalc(1)*1e6, zPtsCalc(end)*1e6]);
ylabel ('Numeric Error[%]','FontSize',axisFont)

end

