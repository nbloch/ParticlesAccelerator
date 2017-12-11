function [ fig ] = plotPhaseSpace( focusedMask, hitMask, numOfParticles, exitR, entryRvec, ...
                                    startGamma, startBetaR, endGamma, endBetaR, Y, notFocusedMask, plotNotFocused)
%PLOTPHASESPACE Summary of this function goes here
%   Detailed explanation goes here
passMask = logical(focusedMask + notFocusedMask);
entryPhase = startGamma.*startBetaR;
endPhase   = endGamma.*endBetaR;
if(sum(focusedMask) == 0 )
    fig = 0;
    return
end

fig = figure();
    legPlot(1) = plot(NaN,NaN,'r+');    hold on;
    legPlot(2) = plot(NaN,NaN,'b+');

    
    plot(exitR(focusedMask),  endPhase(focusedMask), 'r+')
    plot(entryRvec(focusedMask),entryPhase(focusedMask), 'b+')
 
    i=1;
    if(plotNotFocused)
        legPlot(3) = plot(NaN,NaN, '-');
        legPlot(4) = plot(NaN,NaN,'--');
        plot(exitR(notFocusedMask), endPhase(notFocusedMask), 'r+')
        plot(entryRvec(notFocusedMask), entryPhase(notFocusedMask), 'b+')
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
        xl = max(abs([min(min(exitR(passMask)),     min(entryRvec (passMask)))  max(max(exitR(passMask)),     max(entryRvec(passMask)))]));
        yl = max(abs([min(min(min(endPhase(passMask))), min(entryPhase(passMask)))  max(max(max(endPhase(passMask))), max(entryPhase(passMask)))]));
        xlim(1.1.*[-xl xl]);
        ylim(1.1.*[-yl yl]); 
    else
        while i <= numOfParticles
                while ((focusedMask(i) == 0) && i<= numOfParticles)
                    i=i+1;
                    if(i>numOfParticles); break; end
                end
                if (i>numOfParticles); break; end
                plot([exitR(i), entryRvec(i)],[endPhase(i), entryPhase(i)])
                i=i+1;
        end
        hold off
        leg = legend (legPlot,'Exit Phase', 'Start Phase');
        leg.Location = 'northeast';
        xl = max(abs([min(min(exitR(focusedMask)),     min(entryRvec (focusedMask)))  max(max(exitR(focusedMask)),     max(entryRvec(focusedMask)))]));
        yl = max(abs([min(min(min(endPhase(focusedMask))), min(entryPhase(focusedMask)))  max(max(max(endPhase(focusedMask))), max(entryPhase(focusedMask)))]));
        xlim(1.1.*[-xl xl]);
        ylim(1.1.*[-yl yl]);   
    end 
    xlabel('Entry/Exit R')
    ylabel('\gamma\beta_r')
    title('Particle Phase Space')

end

