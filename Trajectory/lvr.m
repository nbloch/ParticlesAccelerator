close all;
clear all;
clc;

c0 = 3e8;
e0 = -1.60217662e-19;
eM = 9.10938356e-31;          

%device Parameters
deviceRadius              = 5e-3;
lensPreOffset             = 2e-3;
lensPostOffset            = 2e-3;
distanceBetweenElectrodes = 1e-3;
leftElectrodeRadius       = 1e-3;
rightElectrodeRadius      = 1e-3;
electrodeWidth            = 0.01e-6;
repetitions               = 1;
VaLeft                    = -15e3;
VaRight                   = 15e3;
globalVa                  = 15e3;

%Resolution Parameters     
M                         = 730;
N                         = 290;     
rPts                      = 1000;
zPts                      = 1000;


%Particle Parameters
electrodeProximityThresh  = 0.01e-6;
q                         = 1.60217662e-19;
m                         = 9.10938356e-31;
simulateElectron          = true;

%Multi Particle
numOfParticles    = 1000;
beamInitialRadius = 0.5e-6;
maxInitMoment     = 1e-3;
Ek                = 10*1e3;

%Single Particle - irrelevant for this sim
axialEntryVel  = 0.2*c0;
radialEntryVel = 0.05*c0;
entryR         = -0.5e-3;
exitRthresh    = 0.25e-6;

deviceLength = repetitions*distanceBetweenElectrodes + lensPreOffset + lensPostOffset;

% =========================
% Build Lens & Potential calc
% =========================
[V, zGrid, rGrid, MSE, Rq, Zq, Rb, Zb, MLeft, MRight, NLeft, Qvec] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                              deviceRadius, distanceBetweenElectrodes, N, M, repetitions, rPts, zPts, lensPreOffset, lensPostOffset);


electordesZ = unique(Zq);
entryZ   = zGrid(1,1);

% =========================
% Low Relativistic Particle
% =========================
[ Zl, Xl, Vzl, Vxl, ~] =  relativeParticleTrajectory(V, zGrid, rGrid, q, m, axialEntryVel,...
                                 radialEntryVel, entryR, entryZ, electrodeProximityThresh,...
                                 deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);
                                                
[ Zr, Xr, Vzr, Vxr, ~ ] = ParticleTrajectory(V, zGrid, rGrid, q, m, axialEntryVel,...
                                 radialEntryVel, entryR, entryZ, electrodeProximityThresh,...
                                 deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);

% =========================
% High Relativistic Particle
% =========================

axialEntryVel       = 0.9  * c0;
radialEntryVel      = 0.1  * c0;                                           
[ Zlr, Xlr, Vzlr, Vxlr, ~] =  relativeParticleTrajectory(V, zGrid, rGrid, q, m, axialEntryVel,...
                                 radialEntryVel, entryR, entryZ, electrodeProximityThresh,...
                                 deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);

[ Zrr, Xrr, Vzrr, Vxrr, ~ ] = ParticleTrajectory(V, zGrid, rGrid, q, m, axialEntryVel,...
                                             radialEntryVel, entryR, entryZ, electrodeProximityThresh,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);

%TODO: add logic to calculate without the complete code
scale = 1e3;

figure()
subplot (1,2,1)
plot(Zr(1:10:end)*scale,Xr(1:10:end)*scale,'-o',Zl(1:10:end)*scale,Xl(1:10:end)*scale,'-+')
title('Energy Consercvation Approach (by Levi) vs. Force Analysis (by Reiser) - Not Relativistic Particle')
xlabel('Z[mm]')
ylabel('R[mm]')
legend('Force(Reiser)', 'Energy(Levi)')
dist = zeros(1,length(Xr));
for i=1:length(Xr)
    [M, Ix] = min(abs(Zl - Zr(i)));
    dist(i) = sqrt( (Zl(Ix)-Zr(i))^2 + (Xl(Ix)-Xr(i))^2);
end
subplot (1,2,2)
plot(Zr*scale, dist*scale)
title('Distance Between 2 Methods Trajectories');
xlabel('Z [mm]')
ylabel('Distance [mm]')

figure()
subplot (1,2,1)
plot(Zrr(1:10:end)*scale,Xrr(1:10:end)*scale,'-o',Zlr(1:10:end)*scale,Xlr(1:10:end)*scale,'-+')
title('Energy Conservation Approach (by Levi) vs. Force Analysis (by Reiser) - Relativistic Particle')
xlabel('Z[mm]')
ylabel('R[mm]')
legend('Force(Reiser)', 'Energy(Levi)')
dist = zeros(1,length(Zrr));
for i = 1:length(Zrr)
    [M, Ix] = min(abs(Zlr - Zrr(i)));
    dist(i) = sqrt( (Zrr(i)-Zlr(Ix))^2 + (Xrr(i)-Xlr(Ix))^2);
end
subplot (1,2,2)
plot(Zrr*scale, dist*scale)
title('Distance Between 2 Methods Trajectories');
xlabel('Z [mm]')
ylabel('Distance [mm]')


figure()
contourf(zGrid*1e3, rGrid*1e3, V/(1e3), 30);
title('Particle Trajectory In Electrostatic Lens - Relativistic');
xlabel('Z axis [mm]')
ylabel('R axis [mm]')
hold on
shading interp
h = colorbar;
h.Label.String = 'V_a [ KV ]';
caxis([min(VaLeft, VaRight) max(VaLeft, VaRight)]/1e3);
colormap(parula);
plot(Zrr*1e3,Xrr*1e3,'r')

figure()
contourf(zGrid*1e3, rGrid*1e3, V/(1e3), 30);
title('Particle Trajectory In Electrostatic Lens - Not Relativistic');
xlabel('Z axis [mm]')
ylabel('R axis [mm]')
hold on
shading interp
h = colorbar;
h.Label.String = 'V_a [ KV ]';
caxis([min(VaLeft, VaRight) max(VaLeft, VaRight)]/1e3);
colormap(parula);
plot(Zr*1e3,Xr*1e3,'r')