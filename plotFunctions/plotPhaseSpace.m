function [ fig ] = plotPhaseSpace( focusedMask, hitMask, numOfParticles, exitR, entryRvec, ...
                                    startGamma, startBetaR, endGamma, endBetaR, Y)
%PLOTPHASESPACE Summary of this function goes here
%   Detailed explanation goes here
passMask = ~hitMask;
notFocusedMask = logical((passMask)-focusedMask);
entryPhase = startGamma.*startBetaR;
endPhase   = endGamma.*endBetaR;

fig = figure();
    legPlot(1) = plot(NaN,NaN,'r+');    hold on;
    legPlot(2) = plot(NaN,NaN,'b+');
    legPlot(3) = plot(NaN,NaN, '-');
    legPlot(4) = plot(NaN,NaN,'--');
    
    plot(exitR(focusedMask),  endPhase(focusedMask), 'r+')
    plot(entryRvec(focusedMask),entryPhase(focusedMask), 'b+')
    plot(exitR(notFocusedMask), endPhase(notFocusedMask), 'r+')
    plot(entryRvec(notFocusedMask), entryPhase(notFocusedMask), 'b+')
    
 
    i=1;
    while i <= numOfParticles
            while ((focusedMask(i) == 0) && (notFocusedMask(i) == 0) && i<= numOfParticles)
                i=i+1;
                if(i>numOfParticles); break; end
            end
            if (i>numOfParticles); break; end
            if (focusedMask(i) == 1)
                plot([exitR(i), entryRvec(i)],[endPhase(i), entryPhase(i)])
            else
                plot([exitR(i), entryRvec(i)],[endPhase(i), entryPhase(i)],'--')
            end
            i=i+1;
    end
    hold off
    leg = legend (legPlot,'Exit Phase', 'Start Phase', 'Focused', 'Not Focused');
    leg.Location = 'northeast';
    xl = max(abs([min(min(min(Y(passMask,:))),     min(entryRvec (passMask)))  max(max(max(Y(passMask,:))),     max(entryRvec(passMask)))]));
    yl = max(abs([min(min(min(endPhase(passMask))), min(entryPhase(passMask)))  max(max(max(endPhase(passMask))), max(entryPhase(passMask)))]));
    xlim(1.1.*[-xl xl]);
    ylim(1.1.*[-yl yl]);     
    xlabel('Entry/Exit R')
    ylabel('\gamma\beta_r')
    title('Particle Phase Space')

end

