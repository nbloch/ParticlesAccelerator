%Generated just to plot problematic particles who are going backwards
%during the simulation
if isdir('./problematicParticles/')
    rmdir('./problematicParticles','s');
end
mkdir ('./problematicParticles/')

k=1;
j=1;

for i = 1:numOfParticles
   if (min(U(i,:)) < 0) %if particle started a motion backwards 
        if (mod(k,10) == 1) 
            k=1;
            probFig = figure();
            j=j+1;

            ax1 = axes('Position',[0 0.05 0.5 0.815],'Visible','off');
            ax2 = axes('Position',[0.17 0.1 0.8 0.8],'Visible','off');
            axes(ax2)
            hq = quiver(zGrid, rGrid, q*Ez,q*Er, 'MaxHeadSize', 100,'AutoScale', 'on', 'AutoScaleFactor', 10, 'LineWidth', 1);           %get the handle of quiver
            hold on;
            title (sprintf('Problematic Particle: %d - Trajectory Over Force Applied On A Single Electron( F=q*E)', i));
            xlabel('Z axis [mm]');
            ylabel('R axis [mm]');
            xlim([zGrid(1,1) zGrid(1,end)]);
           
            axes(ax1)

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
            };
            text(.025, 0.55, params_str);
            axes(ax2)
            legstr{1} = 'Force applied on a single electron'; 
 
        end
        
        if (simulatePhaseSpace)
            plot(X(i,1:trajectoryLen(i)),Y(i,1:trajectoryLen(i)));
            legstr{k+1} = sprintf('particle #: %d, V_z_-_i_n = %.2f[c],\nV_r_-_i_n = %.2f[c], R_i_n = %.2f[mm]',...
                                   i ,axialEntryVelVec(i)/c0, radialEntryVelVec(i)/c0, entryRvec(i)*1e3);
        else 
            plot(X,Y);
            legstr{2}   = sprintf('particle #: %d, V_z_-_i_n = %.2f[c],\nV_r_-_i_n = %.2f[c], R_i_n = %.2f[mm]',...
                                   1 ,axialEntryVelocity/c0, radialEntryVelocity/c0, entryR*1e3);
            k=0;
        end
        
        if (mod(k,10) == 0)
            legend(legstr, 'Location', 'eastoutside');
            hold off;
            savefig(probFig, sprintf('./problematicParticles/problematicParticles-%d', i))
            close(probFig)
        end
        k = k+1;
   end    
end

if (k~=1)
    legend(legstr, 'Location', 'eastoutside');
    savefig(probFig,sprintf('./problematicParticles/problematicParticles-%d', numOfParticles))
end
    
