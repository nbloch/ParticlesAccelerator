function [] = phaseSpaceVideo( Z, X, Vz, Vx, entryRvec, zGrid, startBetaR, startGamma, params,lambda0)
%PHASESPACEVIDEO Summary of this function goes here
%   Detailed explanation goes here
c0 = 3e8;

gamma = 1./sqrt(1-((Vx.^2+Vz.^2)./(c0^2)));
phase = gamma.*Vx./c0;
startPhase = startGamma.*startBetaR;

videoFWriter = vision.VideoFileWriter('phaseSpace.avi','FrameRate',10,'FileFormat','AVI');
vidParams = {' '; '------Shifting-------'; ' ';' ';' ';};
vidParams = [params; vidParams];   
xl = max(abs([min([min(min(X)), min(entryRvec)]), max([max(max(X)), max(entryRvec)])]));
yl = max(abs([min([min(min(phase)), min(startPhase)])  max([max(max(phase)), max(startPhase)])]));
xlimits = 1.1.*[-xl xl];
ylimits = 1.1.*[-yl yl];

rows = 1:size(X,1); 
for j=1:5:length(zGrid(1,:))
    phaseSpaceVidFig = figure();
    set(gcf, 'Position', [0 0 1920 1080])
    ax1 = axes('Position',[0 0.05 0.5 0.815],'Visible','off');
    ax2 = axes('Position',[0.17 0.1 0.8 0.8],'Visible','off');
    axes(ax2)
    
    %finding the points and calculating the emittance
    z = zGrid(1,j);
    [~, cols]= min(abs(Z-z),[],2);
    idxs = sub2ind(size(X), rows, cols');
    emittance = getEmittance( X(idxs), phase(idxs), gamma(idxs), lambda0 );
    vidParams{end-2} = ['z: ',num2str(z)];
    vidParams{end-1} = ['\epsilon_c: ',num2str(emittance)];
    vidParams{end}   = ['P_x(max): ',num2str(max(phase(idxs)))];
    
    %plotting
    plot(entryRvec(:), startPhase(:), 'b+');
    hold on;
    plot(X(idxs), phase(idxs), 'r+')
    for i =1:rows(end)
        plot([X(idxs(i)), entryRvec(i)],[phase(idxs(i)), startPhase(i)])
    end
    title('Changing Phase Space as Beam Spread')
    legend ({'Exit Phase'; 'Start Phase'}, 'Location', 'northeast')
    xlabel('Entry/Exit R')
    ylabel('\gamma\beta_r')
    xlim(xlimits);
    ylim(ylimits);
    axes(ax1);
    text(.025,0.55, vidParams);
    
    %record the fig
    F = getframe(phaseSpaceVidFig);
    step(videoFWriter,F.cdata);
    close(phaseSpaceVidFig);
end
fprintf('\nDone Video\n')
release(videoFWriter);
end