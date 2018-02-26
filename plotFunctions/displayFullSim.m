function [fig] = displayFullSim(VaLeft, VaRight, zGrid, rGrid, simulatePhaseSpace, numOfParticles, ...
                                trajectoryLen, V, X, Y, simulateTrajectory, Mbleft, M, repetitions, ...
                                Zb, Rb, params_str)
% Display Lens and Trajectories

dimFactor = 1e6;

fig = figure();
ax1 = axes('Position',[0 0.05 0.5 0.815],'Visible','off');
ax2 = axes('Position',[0.17 0.1 0.8 0.8],'Visible','off');
axes(ax2)
contourf(zGrid*dimFactor, rGrid*dimFactor, V, 30);
title('Particle Trajectory In Electrostatic Lens');
ax2.TitleFontSizeMultiplier = 2;
xlabel('Z axis [\mum]')
ylabel('X axis [\mum]')
hold on
shading interp
h = colorbar;
h.Label.String = 'V_a [ V ]';
caxis([min(VaLeft, VaRight) max(VaLeft, VaRight)]);
colormap(parula);

if (simulatePhaseSpace)
     for i =1:numOfParticles
         plot(X(i,1:trajectoryLen(i))*dimFactor,Y(i,1:trajectoryLen(i))*dimFactor)
     end
elseif (simulateTrajectory)
     plot(X*dimFactor,Y*dimFactor,'r')
end
 
 %We add the electrodes shape
 if VaLeft<VaRight
     leftColor = [0 0 0];
     rightColor = [1 0 0];
 else 
     leftColor = [1 0 0];
     rightColor = [0 0 0];
 end
 leftIdx = [1:Mbleft];
 rightIdx = [Mbleft+1:M];
 for i = 2:repetitions
     leftIdx = [leftIdx, (i-1)*M+1:(i-1)*M+Mbleft];
 end
 for i = 1:repetitions
     rightIdx = [rightIdx, (i-1)*M+Mbleft+1:i*M];
 end
 leftIdx = [leftIdx, repetitions*M+1:repetitions*M+Mbleft];
 
%left electrode
scatter([Zb(leftIdx)*dimFactor; Zb(leftIdx)*dimFactor],[Rb(leftIdx)*dimFactor; -Rb(leftIdx)*dimFactor], 11, leftColor, 'square', 'filled');
%right electrode
scatter([Zb(rightIdx)*dimFactor; Zb(rightIdx)*dimFactor],[Rb(rightIdx)*dimFactor; -Rb(rightIdx)*dimFactor], 11, rightColor, 'square', 'filled');

hold off
if (repetitions > 2)
    axis equal
end
axes(ax1);
text(.025,0.55, params_str,'FontSize',11);
axes(ax2)

end

