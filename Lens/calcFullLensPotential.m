function [V, zGrid, rGrid, MSE, fig, Rq, Zq, Rb, Zb, MLeft, MRight] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                              deviceRadius, distanceBetweenElectrodes, N, M, use_bessel, repetitions, deviceLength,...
                                              rPts, zPts)
%-----------------------------
%Calculate Charge Distribution
%-----------------------------

[Q, Rq, Zq, Rb, Zb, MSE, NLeft, NRight, lastZoff, MLeft, MRight] = getChargesRepetitive(VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                              deviceRadius, distanceBetweenElectrodes, N, M, use_bessel, repetitions );
fprintf('Done Getting Charges\n');

%---------------------
%Plotting
%---------------------
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
    hold off;
    title('Left Electrode Charge Distribution');
    xlabel('R [m]');
    ylabel('C [Farad]');
    legend(legstr,'Location', 'best');
    ylim([minlim maxlim]);
    subplot(1,2,2)
    for i = 1:repetitions
        plot(Rq((i-1)*N+1+NLeft:i*N), Q((i-1)*N+1+NLeft:i*N)/(sign(VaRight)*Vdiff));
        hold on;
    end
    hold off;
    title('Right Electrode Charge Distribution');
    xlabel('R [m]');
    ylabel('C [Farad]');
    legend(legstr(1:end-1),'Location', 'best');
    ylim([minlim maxlim]);

%---------------------
%Create Grid
%---------------------

r = linspace(0, deviceRadius, rPts);
z = linspace (0, deviceLength, zPts);
z = [-flip(z), z];
z(zPts) = [];

%{
if repetitions > 1
    for i = 1:repetitions-1
        tmp = linspace (0, deviceLength, zPts);
        tmp = [-flip(tmp), tmp];
        tmp(zPts) = [];
        z = [z-deviceLength , tmp+deviceLength];
    end
end
%}

[z_mat, r_mat] = meshgrid(z,r);

flippedR = -flip(r_mat);
flippedR(1,:) = [];
rGrid = [flippedR ; r_mat];
zGrid = [z_mat; z_mat];
zGrid(end,:) = [];

%----------------------------------
%Calculating field Of the Full Lens
%----------------------------------
r_mat = repmat(r_mat,1,1,length(Rq));
z_mat = repmat(z_mat,1,1,length(Rq));
R = repmat(reshape(Rq,1,1,length(Rq)),size(r_mat,1), size(r_mat,2));
Z = repmat(reshape(Zq,1,1,length(Zq)),size(r_mat,1), size(r_mat,2));
Q  =repmat(reshape(Q,1,1,length(Q)),size(r_mat,1), size(r_mat,2));


if (use_bessel)
    V = sum(Q.*BesselRingPotential(r_mat, z_mat, R, Z),3);
else
   V = sum(Q.*spatial_green(r_mat, z_mat, R, Z),3);
end
fprintf('Done Computing Full Space Potential\n');

Vflipped = flip(V);
Vflipped(1,:) = [];
V = [Vflipped; V];

%We want to avoid very high value points that may make the colormap less readable
% Vmax = max(VaLeft, VaRight);
% Vmin = min(VaLeft, VaRight);
% V(V > Vmax) = Vmax;
%We want to avoid very high value points that may make the colormap less readable
% Vmax = max(VaLeft, VaRight);
% Vmin = min(VaLeft, VaRight);
% V(V > Vmax) = Vmax;
% V(V<Vmin) = Vmin;

end
