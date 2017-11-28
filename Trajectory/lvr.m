VaLeft                      = -15e3;
VaRight                     =  15e3;
electrodeWidth              =  10e-6;
leftElectrodeRadius         =  1e-3;
rightElectrodeRadius        =  1e-3;     
deviceRadius                =  5e-3;
distanceBetweenElectrodes   =  1e-3;
N                           =  250;
M                           =  600;
use_bessel                  =  0;
rPts                        =  1000;
zPts                        =  1000; 
sideOffset                  =  2e-3; 
repetitions  = 1;

c0 = 3e8;


simulatePhaseSpace       = 0;
useRelativeTrajectory    = 1;
q                        = -1.60217662e-19;
m                        = 9.10938356e-31;
axialEntryVelocity       = 0.1  * c0;
radialEntryVelocity      = 0.05 * c0;

entryR                   = 0; 
electrodeProximityThresh = 10e-6;
exitRthresh              = 0.25*1e-3;
lowAxialVel              = 0.1*c0;
highAxialVel             = 0.1*c0;
lowRadialVel             = (1e-6)*c0;
highRadialVel            = (1e-6)*c0;
numOfParticles           = 1;



deviceLength = repetitions*distanceBetweenElectrodes + sideOffset; 

[V, zGrid, rGrid, MSE, fig(2), Rq, Zq,  Rb, Zb ] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                                    deviceRadius, distanceBetweenElectrodes, N, M, use_bessel, repetitions, deviceLength,...
                                                    rPts, zPts);

[ XLR, YLR] = calcTrajectory( simulatePhaseSpace,useRelativeTrajectory, V, zGrid, rGrid, q, m, axialEntryVelocity,...
                                             radialEntryVelocity, entryR, Zq ,...
                                             electrodeProximityThresh, exitRthresh,...
                                             lowAxialVel, highAxialVel, lowRadialVel, highRadialVel, numOfParticles,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius);
                                                

useRelativeTrajectory    = 0;
[ XR, YR] = calcTrajectory( simulatePhaseSpace,useRelativeTrajectory, V, zGrid, rGrid, q, m, axialEntryVelocity,...
                                             radialEntryVelocity, entryR, Zq ,...
                                             electrodeProximityThresh, exitRthresh,...
                                             lowAxialVel, highAxialVel, lowRadialVel, highRadialVel, numOfParticles,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius);
                                         
axialEntryVelocity       = 0.9  * c0;
radialEntryVelocity      = 0.1  * c0;                                         
useRelativeTrajectory    = 1;                                         
[ XL, YL] = calcTrajectory( simulatePhaseSpace,useRelativeTrajectory, V, zGrid, rGrid, q, m, axialEntryVelocity,...
                                             radialEntryVelocity, entryR, Zq ,...
                                             electrodeProximityThresh, exitRthresh,...
                                             lowAxialVel, highAxialVel, lowRadialVel, highRadialVel, numOfParticles,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius);
                                                

useRelativeTrajectory    = 0;
[ X, Y] = calcTrajectory( simulatePhaseSpace,useRelativeTrajectory, V, zGrid, rGrid, q, m, axialEntryVelocity,...
                                             radialEntryVelocity, entryR, Zq ,...
                                             electrodeProximityThresh, exitRthresh,...
                                             lowAxialVel, highAxialVel, lowRadialVel, highRadialVel, numOfParticles,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius);                                         


%TODO: add logic to calculate without the complete code


figure()
subplot (1,2,1)
plot(X(1:10:end),Y(1:10:end),'-o',XL(1:10:end),YL(1:10:end),'-+')
title('Levi eq vs. Reiser Eq - Not Relativistic Particle')
xlabel('Z[mm]')
ylabel('R[mm]')
legend('Reiser', 'Levi')
dist = zeros(1,length(X));
for i=1:length(X)
    [M, Ix] = min(abs(XL - X(i)));
    dist(i) = sqrt( (X(i)-XL(Ix))^2 + (Y(i)-YL(Ix))^2);
end
subplot (1,2,2)
plot(X, dist)
title('Distance Between 2 Methods Trajectories');
xlabel('Z [mm]')
ylabel('Distance [m]')

figure()
subplot (1,2,1)
plot(XR(1:10:end),YR(1:10:end),'-o',XLR(1:10:end),YLR(1:10:end),'-+')
title('Levi eq vs. Reiser Eq - Relativistic Particle')
xlabel('Z[mm]')
ylabel('R[mm]')
legend('Reiser', 'Levi')
dist = zeros(1,length(X));
for i=1:length(XR)
    [M, Ix] = min(abs(XLR - XR(i)));
    dist(i) = sqrt( (XR(i)-XLR(Ix))^2 + (YR(i)-YLR(Ix))^2);
end
subplot (1,2,2)
plot(XR, dist)
title('Distance Between 2 Methods Trajectories');
xlabel('Z [mm]')
ylabel('Distance [m]')


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
plot(XR*1e3,YR*1e3,'r')

figure()
contourf(zGrid*1e3, rGrid*1e3, V/(1e3), 30);
title('Particle Trajectory In Electrostatic Lens - Nor Relativistic');
xlabel('Z axis [mm]')
ylabel('R axis [mm]')
hold on
shading interp
h = colorbar;
h.Label.String = 'V_a [ KV ]';
caxis([min(VaLeft, VaRight) max(VaLeft, VaRight)]/1e3);
colormap(parula);
plot(X*1e3,Y*1e3,'r')