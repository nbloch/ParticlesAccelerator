function [V, zGrid, rGrid, MSE, Rq, Zq, Rb, Zb, MLeft, MRight, NLeft, Qvec] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                              deviceRadius, distanceBetweenElectrodes, N, M, repetitions, rPts, zPts, lensPreOffset, lensPostOffset)
%-----------------------------
%Calculate Charge Distribution
%-----------------------------

[Q, Rq, Zq, Rb, Zb, MSE, NLeft, NRight, lastZoff, MLeft, MRight] = getChargesRepetitive(VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                              deviceRadius, distanceBetweenElectrodes, N, M, repetitions );
fprintf('Done Getting Charges\n');
Qvec = Q;

%---------------------
%Create Grid
%---------------------
leftLength = -lensPreOffset-(repetitions*distanceBetweenElectrodes)/2;
rightLength = lensPostOffset+(repetitions*distanceBetweenElectrodes)/2;

r = linspace(0, deviceRadius, rPts);
z = linspace (leftLength, rightLength, 2*zPts);
% z = [-flip(z), z];
% z(zPts) = [];

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

V = sum(Q.*spatial_green(r_mat, z_mat, R, Z),3);

fprintf('Done Computing Full Space Potential\n');

Vflipped = flip(V);
Vflipped(1,:) = [];
V = [Vflipped; V];
end
