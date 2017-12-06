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
% fig = 0;
% if (plotResults)
    fig = displayChargeDistribution( repetitions, Rq, Q, N, NLeft, VaLeft, VaRight );
% end 

%---------------------
%Create Grid
%---------------------

r = linspace(0, deviceRadius, rPts);
z = linspace (0, deviceLength, zPts);
z = [-flip(z), z];
z(zPts) = [];

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
