
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
rPts                        =  200;
zPts                        =  200; 
sideOffset                  =  2e-3; 

MSE = zeros(3,7);
for i = 1:7
    deviceLength = i*distanceBetweenElectrodes + sideOffset; 
    [~, ~, ~, MSE(1,i)] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                                        deviceRadius, distanceBetweenElectrodes, N, M, use_bessel, i, deviceLength,...
                                                        rPts, zPts);
end 


N                           =  300;
M                           =  1000;
for i = 1:7
    deviceLength = i*distanceBetweenElectrodes + sideOffset; 
    [~, ~, ~, MSE(2,i)] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                                        deviceRadius, distanceBetweenElectrodes, N, M, use_bessel, i, deviceLength,...
                                                        rPts, zPts);
end 

N                           =  400;
M                           =  1600;
for i = 1:7
    deviceLength = i*distanceBetweenElectrodes + sideOffset; 
    [~, ~, ~, MSE(3,i)] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                                        deviceRadius, distanceBetweenElectrodes, N, M, use_bessel, i, deviceLength,...
                                                        rPts, zPts);
end 

                                        
figure()
for i = 1:3
    plot(1:7, MSE(i, :), '-o')
    hold on
end
title('MSE vs, Number Of Unit Cells')
xlabel('Number of Unit Cells')
ylabel('MSE')
hleg = legend('M = 600, N=250', 'M = 1000, N=300', 'M = 1600, N=400', 'Position', 'northeast');
htitle = get(hleg, 'Title');
set(htitle,'String', 'Number Of Variables Per Cell M- Boundaries, N- Charges');




                                                
                                                