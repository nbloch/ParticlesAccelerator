function [ V, X, Y, U, W, MSE, fig, zGrid, rGrid, Rq, Zq, Rb, Zb, Ez, Er, Mbleft, Mbright, focused_particles ] =...
    FullEinzelSim(simPars)
%NOTE: The parameters "calcPotential, V, zGrid, rGrid, Rq, Zq, MSE" are
% used for device potential caching. if calcPotential is set to true, the
% other parameters are unused and can be left empty

%Load struct parameters into workspace
save('params-tmp.mat', '-struct', 'simPars');
load('params-tmp.mat');
delete('params-tmp.mat');
clear params;
%END

close all;
fig = gobjects(1,4);
addpath(genpath([pwd]));
%--------------------------------
%Params For trajectory
%--------------------------------

deviceLength = repetitions*distanceBetweenElectrodes + sideOffset;

if(simulateElectron)
   q =  -1.60217662e-19;
   m = 9.10938356e-31;
end
c0 = 3e8;

if (~simulatePhaseSpace)
    numOfParticles = 1;
end

% if(useAngle)
%     %For one particle
%     axialEntryVelocity  = entryVel*cos(deg2rad(entryAngle));
%     radialEntryVelocity = entryVel*sin(deg2rad(entryAngle));
%     
% %     %For many particles
% %     lowAxialVel   = entryVel*cos(deg2rad(entryAngle));
% %     highAxialVel  = entryVel;
% %     lowRadialVel  = (1e-6)*c0;
% %     highRadialVel = entryVel*sin(deg2rad(entryAngle));
% end


%--------------------------------
%Calculating Operating Point
%--------------------------------
if(dimensionsOptimization)
    tic;
    [ N, M, ~, fig(1) ] = FindOperatingPoint(VaLeft, VaRight, electrodeWidth,leftElectrodeRadius,...
                                         rightElectrodeRadius, deviceRadius, distanceBetweenElectrodes,...
                                         use_bessel, repetitions, deviceLength,  convergeTh, growthTh );
    toc;
end


%--------------------------------
%Calculating Potential
%--------------------------------
if(calcPotential)
      tic;
      [V, zGrid, rGrid, MSE, fig(2), Rq, Zq,  Rb, Zb, Mbleft, Mbright ] = calcFullLensPotential( VaLeft, VaRight, electrodeWidth, ...
                                                    leftElectrodeRadius, rightElectrodeRadius, deviceRadius, ...
                                                    distanceBetweenElectrodes, N, M, use_bessel, repetitions, deviceLength,...
                                                    rPts, zPts);
      toc;
end

%--------------------------------
%Calculating Particle Trajectory
%--------------------------------
Rq = [Rq, -Rq];
Zq = [Zq,  Zq];
X = []; Y = []; U = []; W = [];
z_step = zGrid(1,2) - zGrid(1,1);
r_step = rGrid(2,1) - rGrid(1,1);
[Ez, Er] = gradient(-V, z_step, r_step);
if(simulateTrajectory)
  tic;
  [ X, Y, U, W, fig(3),parameters, trajectoryLen, focused_particles] = calcTrajectory( simulatePhaseSpace,...
                                             useRelativeTrajectory, V, zGrid, rGrid, q, m, entryVel,...
                                             entryAngle, axialEntryVelocity, radialEntryVelocity, entryR, Zq ,...
                                             electrodeProximityThresh, beamInitialRadius,...
                                             lowAxialVel, highAxialVel, lowRadialVel, highRadialVel, numOfParticles,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, VaLeft, VaRight, Ez, Er, ...
                                             recordPhaseSpace, useAngle, beamInitialRadius);
  toc;
end

%--------------------------------
%Plotting
%--------------------------------
params_str = { '--------Device Parameters--------';
              ['Electrode Width: ', num2str(electrodeWidth*1e6),'[\mum]'];
              ['R_a_p_-_l: ', num2str(leftElectrodeRadius*1e3),'[mm]'];
              ['R_a_p_-_r: ', num2str(rightElectrodeRadius*1e3),'[mm]'];
              ['R_L_e_n_s: ', num2str(deviceRadius*1e3),'[mm]'];
              ['d_a_p: ', num2str(distanceBetweenElectrodes*1e3),'[mm]'];
              ['l_L_e_n_s: ', num2str(deviceLength*1e3),'[mm]'];
               ' ';
               '-------Electric Parameters-------';
              ['V_a_-_L_e_f_t: ', num2str(VaLeft),'[V]'];
              ['V_a_-_R_i_g_h_t: ', num2str(VaRight),'[V]'];
               ' ';
               '------Resolution Parameters------';
              ['MSE: ', num2str(MSE),'[%]'];
              ['M_b: ', num2str(M)];
              ['N_c: ', num2str(N)];
		};

if(simulateTrajectory)
    params_str = [params_str ; parameters];
end

fig(4) = displayFullSim(VaLeft, VaRight, zGrid, rGrid, simulatePhaseSpace, numOfParticles, ...
                                trajectoryLen, V, X, Y, simulateTrajectory, Mbleft, M, repetitions, ...
                                Zb, Rb, params_str)
                            
% fig(4) = figure();
% ax1 = axes('Position',[0 0.05 0.5 0.815],'Visible','off');
% ax2 = axes('Position',[0.17 0.1 0.8 0.8],'Visible','off');
% axes(ax2)
% contourf(zGrid*1e3, rGrid*1e3, V/(1e3), 30);
% 
% 
% title('Particle Trajectory In Electrostatic Lens');
% xlabel('Z axis [mm]')
% ylabel('R axis [mm]')
% hold on
% shading interp
% h = colorbar;
% h.Label.String = 'V_a [ KV ]';
% caxis([min(VaLeft, VaRight) max(VaLeft, VaRight)]/1e3);
% colormap(parula);
% 
% if (simulatePhaseSpace)
%      for i =1:numOfParticles
%          plot(X(i,1:trajectoryLen(i))*1e3,Y(i,1:trajectoryLen(i))*1e3)
%      end
% elseif (simulateTrajectory)
%      plot(X*1e3,Y*1e3,'r')
% end
%  
%  %We add the electrodes shape
%  if VaLeft<VaRight
%      leftColor = [0 0 0];
%      rightColor = [1 0 0];
%  else 
%      leftColor = [1 0 0];
%      rightColor = [0 0 0];
%  end
%  leftIdx = [1:Mbleft];
%  rightIdx = [Mbleft+1:M];
%  for i = 2:repetitions
%      leftIdx = [leftIdx, (i-1)*M+1:(i-1)*M+Mbleft];
%  end
%  for i = 1:repetitions
%      rightIdx = [rightIdx, (i-1)*M+Mbleft+1:(i-1)*(M+1)];
%  end
%  leftIdx = [leftIdx, repetitions*M+1:repetitions*M+Mbleft];
%  
% %left electrode
% scatter([Zb(leftIdx)*1e3; Zb(leftIdx)*1e3],[Rb(leftIdx)*1e3; -Rb(leftIdx)*1e3], 11, leftColor, 'square', 'filled');
% %right electrode
% scatter([Zb(rightIdx)*1e3; Zb(rightIdx)*1e3],[Rb(rightIdx)*1e3; -Rb(rightIdx)*1e3], 11, rightColor, 'square', 'filled');
% 
% hold off
% if (repetitions > 2)
%     axis equal
% end
% axes(ax1);
% text(.025,0.55, params_str);
% axes(ax2)
end
