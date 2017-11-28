function [ Z, R, U, V, hitElectrode ] = relativeParticleTrajectory(Voltage, z_grid, r_grid, q, m, entryAxialVel,...
                                             entryRadialVel, entryR, entryZ, electrodeProximityThresh,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ)
%PARTICLETRAJECTORY Summary of this function goes here
%   Detailed explanation goes here
   z_step = z_grid(1,2) - z_grid(1,1);
   r_step = r_grid(2,1) - r_grid(1,1);
  [Ez, Er] = gradient(-Voltage, z_step, r_step);
   c0 = 3e8;
% =======================================================================
% Particle Motion Simulate
% =======================================================================
   i  = 1;
%    dt = sqrt(z_step^2 + r_step^2)/(sqrt(entryAxialVel^2+entryRadialVel^2));
   dt = sqrt(z_step^2 + r_step^2)/((sqrt(entryAxialVel^2+entryRadialVel^2))*4);

   q = -q;     %FIXME: Should q be positive or negative?

   g = @(a) 1/(c0^2-(a)^2);
   h = @(u, v) (q/(m*c0))*(c0^2-u^2-v^2)^(3/2);
   fz = @(z,r, u, v) h(u, v)*Ez(getIndexInGrid(r_grid,r),getIndexInGrid(z_grid, z));
   fr = @(z,r, u, v) h(u, v)*Er(getIndexInGrid(r_grid,r),getIndexInGrid(z_grid, z));
   F1 = @(z,r,u,v) g(v)*(g(u)*u*v*fr(z,r,u,v)-fz(z,r,u,v)) ...
      /(1-g(u)*g(v)*(u*v)^2);
   F2 = @(z,r,u,v) g(u)*(g(v)*u*v*fz(z,r,u,v)-fr(z,r,u,v)) ...
      /(1-g(u)*g(v)*(u*v)^2);
   F3 = @(u) u;
   F4 = @(v) v;

   %start Conditions: U=Z', V=R'
   Z(i)     = entryZ;
   R(i)     = entryR;
   U(i)     = entryAxialVel;
   V(i)     = entryRadialVel;


% figure();
% title('Particle Trajectory In Electrostatic Lens');
% xlabel('Z axis [mm]')
% ylabel('R axis [mm]')
% contourf(z_grid,r_grid,Voltage,30);
% hold on
% shading interp
% h = colorbar;
% h.Label.String = 'V [ V ]';
% colormap(parula);
% hold on
   
   
   
   %Runge-Kutta:
   while (Z(i)>=z_grid(1,1) && (Z(i)<=(max(max(z_grid)))) && (abs(R(i))<= max(max(r_grid))))

       %To make the code writing simpler
       Zi=Z(i);
       Ri=R(i);
       Ui=U(i);
       Vi=V(i);

% For more info about RK4 with two variables:
% https://math.stackexchange.com/questions/721076/     Kr0 = F4(V(i));help-with-using-the-runge-kutta-4th-order-method-on-a-system-of-2-first-order-od

%%K0
       Kz0 =  dt*F3(Ui);
       Ku0 = dt*F1(Zi, Ri, Ui, Vi);
       Kr0 = dt*F4(Vi);
       Kv0 = dt*F2(Zi, Ri, Ui, Vi);

%%K1
      Kz1 = dt*F3(Ui+Ku0/2);
      Ku1 = dt*F1(Zi+Kz0/2, Ri+Kr0/2, Ui+Ku0/2, Vi+Kv0/2);
      Kr1 = dt*F4(Vi+Kv0/2);
      Kv1 = dt*F2(Zi+Kz0/2, Ri+Kr0/2, Ui+Ku0/2, Vi+Kv0/2);

%%K2
      Kz2 = dt*F3(Ui+Ku1/2);
      Ku2 = dt*F1(Zi+Kz1/2, Ri+Kr1/2, Ui+Ku1/2, Vi+Kv1/2);
      Kr2 = dt*F4(Vi+Kv1/2);
      Kv2 = dt*F2(Zi+Kz1/2, Ri+Kr1/2, Ui+Ku1/2, Vi+Kv1/2);

%%K3
      Kz3 = dt*F3(Ui+Ku2);
      Ku3 = dt*F1(Zi+Kz2, Ri+Kr2, Ui+Ku2, Vi+Kv2);
      Kr3 = dt*F4(Vi+Kv2);
      Kv3 = dt*F2(Zi+Kz2, Ri+Kr2, Ui+Ku2, Vi+Kv2);

%% i+1 computation
      Z(i+1) = Zi+(Kz0 + 2*Kz1 + 2*Kz2 + Kz3)/6;
      R(i+1) = Ri+(Kr0 + 2*Kr1 + 2*Kr2 + Kr3)/6;
      U(i+1) = Ui+(Ku0 + 2*Ku1 + 2*Ku2 + Ku3)/6;
      V(i+1) = Vi+(Kv0 + 2*Kv1 + 2*Kv2 + Kv3)/6;

       %-----------------------------------
       %check if particle was hit by an electrode.
       %-----------------------------------
       hitElectrode = checkElectrodeProximity( deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ, [R(i) R(i+1)], [Z(i) Z(i+1)], electrodeProximityThresh);
       if (hitElectrode)
           break;
       end
       %-----------------------------------
       i=i+1;
%        plot(Z,R,'-r')
%        pause(0.05)
   end

end
