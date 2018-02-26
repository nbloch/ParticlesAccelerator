function [ fig ] = displayChargeDistribution( repetitions, Rq, Q, N, NLeft, VaLeft, VaRight )
%DISPLAYCHARGEDISTRIBUTION Summary of this function goes here
%   Detailed explanation goes here

for i = 1:repetitions
    legstr{i} = sprintf('Structure = %s',num2str(i));
end
legstr{repetitions+1} = sprintf('Last Electrode');
Vdiff = abs(VaLeft - VaRight);
maxlim = max(abs(Q)/Vdiff);
minlim = min(abs(Q)/Vdiff);

fig = figure();
    subplot(1,2,1)
    
    for i = 1:repetitions+1 %+1 for last electrode
        plot(Rq((i-1)*N+1:(i-1)*N+NLeft), Q((i-1)*N+1:(i-1)*N+NLeft)/(sign(VaLeft)*Vdiff));
        hold on;
    end 
    pbaspect([1 1 1]);
    hold off;
    title('Left Electrode Charge Distribution');
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;
    xlabel('R [m]');
    ylabel('C [Farad]');
    legend(legstr,'Location', 'northwest');
    ylim([minlim maxlim]);
    subplot(1,2,2)
    for i = 1:repetitions
        plot(Rq((i-1)*N+1+NLeft:i*N), Q((i-1)*N+1+NLeft:i*N)/(sign(VaRight)*Vdiff));
        hold on;
    end
    pbaspect([1 1 1]);
    hold off;
    title('Right Electrode Charge Distribution');
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;
    xlabel('R [m]');
    ylabel('C [Farad]');
    legend(legstr(1:end-1),'Location', 'northwest');
    ylim([minlim maxlim]);




end

