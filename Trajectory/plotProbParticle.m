%Generated just to plot problematic particles who are going backwards
%during the simulation

if (simulatePhaseSpace)
    probFig(k) = figure();
else 
    probFig = figure();
    i=1;
end
ax1 = axes('Position',[0 0.05 0.5 0.815],'Visible','off');
ax2 = axes('Position',[0.17 0.1 0.8 0.8],'Visible','off');
axes(ax2)
hq = quiver(zGrid, rGrid, q*Ez,q*Er, 'MaxHeadSize', 100,'AutoScale', 'on', 'AutoScaleFactor', 10, 'LineWidth', 1);           %get the handle of quiver
hold on;
if (simulatePhaseSpace)
    plot(X(i,1:trajectoryLen(i)),Y(i,1:trajectoryLen(i)));
else 
    plot(X,Y);
end
title (sprintf('Problematic Particle: %d - Trajectory Over Force Applied On A Single Electron( F=q*E)', i));
xlabel('Z axis [mm]');
ylabel('R axis [mm]');
xlim([zGrid(1,1) zGrid(1,end)]);
legend('Force applied on a single electron', 'Electron trajectory');
axes(ax1)
if (simulatePhaseSpace)
    params_str = { '--------Device Parameters--------';
  ['R_a_p_-_l: ', num2str(leftElectrodeRadius*1e3),'[mm]'];
  ['R_a_p_-_r: ', num2str(rightElectrodeRadius*1e3),'[mm]'];
  ['R_L_e_n_s: ', num2str(deviceRadius*1e3),'[mm]'];
   ' ';
   '-------Electric Parameters-------';
  ['V_a_-_L_e_f_t: ', num2str(VaLeft),'[V]'];
  ['V_a_-_R_i_g_h_t: ', num2str(VaRight),'[V]'];
   ' ';
   '-------Particle Parameters:------';
  ['q: ', num2str(q/e0),'[e_0]'];
  ['M: ', num2str(m/eM),'[e_M]'];
  ' ';
   '----Particle Entry Parameters----';
  ['V_z_-_i_n: ', num2str(axialEntryVelVec(i)/c0),'[c]'];
  ['V_r_-_i_n: ', num2str(radialEntryVelVec(i)/c0),'[c]'];
  ['R_i_n: ',num2str(entryRvec(i)*1e3),'[mm]'];
  ' ';
   '----Particles Exit Parameters----';
  ['V_z_-_o_u_t: ', num2str(U(i,trajectoryLen(i))/c0),'[c]'];
  ['V_r_-_o_u_t: ', num2str(W(i,trajectoryLen(i))/c0),'[c]'];
  ['R_o_u_t: ',num2str(Y(i,trajectoryLen(i))*1e3),'[mm]'];
  };
else 
    params_str = { '--------Device Parameters--------';
  ['R_a_p_-_l: ', num2str(leftElectrodeRadius*1e3),'[mm]'];
  ['R_a_p_-_r: ', num2str(rightElectrodeRadius*1e3),'[mm]'];
  ['R_L_e_n_s: ', num2str(deviceRadius*1e3),'[mm]'];
   ' ';
   '-------Electric Parameters-------';
  ['V_a_-_L_e_f_t: ', num2str(VaLeft),'[V]'];
  ['V_a_-_R_i_g_h_t: ', num2str(VaRight),'[V]'];
   ' ';
   '-------Particle Parameters:------';
  ['q: ', num2str(q/e0),'[e_0]'];
  ['M: ', num2str(m/eM),'[e_M]'];
  ' ';
   '----Particle Entry Parameters----';
  ['V_z_-_i_n: ', num2str(axialEntryVelocity/c0),'[c]'];
  ['V_r_-_i_n: ', num2str(radialEntryVelocity/c0),'[c]'];
  ['R_i_n: ',num2str(entryR*1e3),'[mm]'];
  ' ';
   '----Particles Exit Parameters----';
  ['V_z_-_o_u_t: ', num2str(X(end)/c0),'[c]'];
  ['V_r_-_o_u_t: ', num2str(W(end)/c0),'[c]'];
  ['R_o_u_t: ',num2str(Y(end)*1e3),'[mm]'];
  ['Hit: ', num2str(hitElectrode)]
  };
end 
text(.025, 0.55, params_str);
axes(ax2)