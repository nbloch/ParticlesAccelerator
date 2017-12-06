function [] = phaseSpaceVideo( q, m, numOfParticles, focused, focusedMask, X, Y, U, W, lowAxialVel, highAxialVel,...
                                lowRadialVel, highRadialVel, entryRvec, zGrid, startBetaR, startGamma)
%PHASESPACEVIDEO Summary of this function goes here
%   Detailed explanation goes here
c0 = 3e8;
e0 =  -1.60217662e-19;
eM = 9.10938356e-31;

gamma = 1./sqrt(1-((W.^2+U.^2)./c0^2));
phase = gamma.*W./c0;
startPhase = startGamma.*startBetaR;

entryR = max(abs(entryRvec))*1e3;
videoFWriter = vision.VideoFileWriter('phaseSpace.avi','FrameRate',10,'FileFormat','AVI');
 vidParams = {' ';
       '------Current Position-------';
       ' ';
       ' ';
       '-----Particle Parameters-----';
      ['q: ', num2str(q/e0),'[e_0]'];
      ['M: ', num2str(m/eM),'[e_M]'];
      ['Particles #: ', num2str(numOfParticles)];
      ['Focused: ', num2str(focused)];
      ' ';
       '--Entry Parameters--';
      ['V_z_-_i_n: [', num2str(lowAxialVel/c0),', ',num2str(highAxialVel/c0),'][c]'];
      ['V_r_-_i_n: [', num2str(lowRadialVel/c0),', ',num2str(highRadialVel/c0),'][c]'];
      ['R_i_n: [', num2str(-entryR),', ',entryR,'][mm]'];
      };
  
xl = max(abs([min(min(min(Y(focusedMask,:))),     min(entryRvec (focusedMask)))  max(max(max(Y(focusedMask,:))),     max(entryRvec(focusedMask)))]));
yl = max(abs([min(min(min(phase(focusedMask,:))), min(startPhase(focusedMask)))  max(max(max(phase(focusedMask,:))), max(startPhase(focusedMask)))]));
xlimits = 1.1.*[-xl xl];
ylimits = 1.1.*[-yl yl];

for j=1:25:length(zGrid(1,:))
    phaseSpaceVidFig = figure();
    ax1 = axes('Position',[0 0.05 0.5 0.815],'Visible','off');
    ax2 = axes('Position',[0.26 0.1 0.7 0.8],'Visible','off');
    axes(ax2)
    z = zGrid(1,j);
    vidParams{3} = sprintf('z: %d',z);
    plot(entryRvec(focusedMask), startPhase(focusedMask), 'b+');
    hold on;
    k=1;
    while k <= numOfParticles
        while ((focusedMask(k) == 0) && k<= numOfParticles)
             k=k+1;
             if(k>numOfParticles); break; end
        end
        if (k>numOfParticles); break; end
        [~, trajIdx] = min(abs(X(k,:)-z));
        plot(Y(k,trajIdx), phase(k,trajIdx), 'r+')
        hold on;
        plot([Y(k,trajIdx), entryRvec(k)],[phase(k,trajIdx), startPhase(k)])
        k=k+1;
    end
    hold off
    title('Changing Phase Space as Beam Spread')
    legend ({'Exit Phase'; 'Start Phase'}, 'Location', 'northeast')
    xlabel('Entry/Exit R')
    ylabel('\gamma\beta_r')
    xlim(xlimits);
    ylim(ylimits);
    axes(ax1);
    text(.025,0.55, vidParams);
    axes(ax2)
    set(gcf, 'Position', [0 0 640, 480])
    F = getframe(phaseSpaceVidFig);
    step(videoFWriter,F.cdata);
    close(phaseSpaceVidFig);
end
fprintf('\nDone Video\n')
release(videoFWriter);

end

