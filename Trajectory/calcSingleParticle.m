function [  Z, X, Vz, Vx, parameters ] = calcSingleParticle( Voltage, zGrid, rGrid, q, m, axialEntryVel, radialEntryVel, entryR, exitRthresh,...
                                             Zq, Ez, Er, VaLeft, VaRight,...
                                             electrodeProximityThresh, deviceRadius, leftElectrodeRadius, rightElectrodeRadius, plotProblematicParticles)
                                         
%------------------------------------%
%-----------Single Particle----------%
%------------------------------------%
tic

c0 = 3e8;
e0 =  -1.60217662e-19;
eM = 9.10938356e-31;
electordesZ = unique(Zq);
entryZ   = zGrid(1,1);

[ Z, X, Vz, Vx, hitElectrode] =  relativeParticleTrajectory(Voltage, zGrid, rGrid, q, m, axialEntryVel,...
                                 radialEntryVel, entryR, entryZ, electrodeProximityThresh,...
                                 deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ);
toc
%------------------------------------%
%---------------Figures--------------%
%------------------------------------%
if (min(Vz) < 0 && plotProblematicParticles) %if particle started a motion backwards          
    displayProbParticle(Z, X, Vz, zGrid, rGrid, Ez, Er, leftElectrodeRadius, rightElectrodeRadius,q, deviceRadius, VaLeft, VaRight,...
                                      m, false, 1, axialEntryVel, radialEntryVel, entryR, 1, Zq);
end
focused = 0;
if (X(end)*1e3 < exitRthresh) && (hitElectrode == 0)
    focused = 1;
end
parameters = {' ';
       '-------Particle Parameters:------';
      ['q: ', num2str(q/e0),'[e_0]'];
      ['M: ', num2str(m/eM),'[e_M]'];
      ' ';
       '----Particle Entry Parameters----';
      ['V_z_-_i_n: ', num2str(axialEntryVel/c0),'[c]'];
      ['V_r_-_i_n: ', num2str(radialEntryVel/c0),'[c]'];
      ['R_i_n: ',num2str(entryR*1e6),'[\mum]'];
      ' ';
       '----Particles Exit Parameters----';
      ['V_z_-_o_u_t: ', num2str(Z(end)/c0),'[c]'];
      ['V_r_-_o_u_t: ', num2str(Vx(end)/c0),'[c]'];
      ['R_o_u_t: ',num2str(X(end)*1e6),'[\mum]'];
      ['Hit: ', num2str(hitElectrode)]
      ['Focused: ', num2str(focused)];
      };

end

