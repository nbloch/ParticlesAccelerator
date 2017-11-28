function [ X, Y, U, W, hitElectrode ] = ParticleTrajectory(Voltage, x_grid, y_grid, q, m, entryAxialVel,...
                                             entryRadialVel, entryR, entryZ, electrodeProximityThresh,...
                                             deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ)
%PARTICLETRAJECTORY Summary of this function goes here
%   Detailed explanation goes here
   x_step = x_grid(1,2) - x_grid(1,1);
   y_step = y_grid(2,1) - y_grid(1,1);
   
  [Ex, Ey] = gradient(-Voltage, x_step, y_step);

% =======================================================================
% Particle Motion Simulate
% =======================================================================
   i  = 1;
   dt = sqrt(x_step^2 + y_step^2)/(sqrt(entryAxialVel^2+entryRadialVel^2));
  
   %start Conditions:
   X(i) = entryZ;
   Y(i) = entryR;
   U(i) = entryAxialVel;
   W(i) = entryRadialVel;
   
   
   %Runge-Kutta:
   while ( X(i)>=x_grid(1,1) && (X(i)<=(max(max(x_grid)))) && (abs(Y(i))<= max(max(y_grid))))

       [~, xCoor] = min(abs(x_grid(1,:)-X(i)));
       [~, yCoor] = min(abs(y_grid(:,1)-Y(i)));
       
       a = calc_RK_order(U(i), W(i),Ex(yCoor, xCoor),Ey(yCoor, xCoor), q, m);
       
       Kx1 = U(i)*dt;
       Ky1 = W(i)*dt;
       Ku1 = a(1)*dt;
       Kw1 = a(2)*dt;

       [~, xCoor] = min(abs(x_grid(1,:)-(X(i)+(Kx1/2))));
       [~, yCoor] = min(abs(y_grid(:,1)-(Y(i)+(Ky1/2))));
       a = calc_RK_order(U(i), W(i),Ex(yCoor, xCoor),Ey(yCoor, xCoor), q, m);

       Kx2 = (U(i)+(Ku1/2))*dt;
       Ky2 = (W(i)+(Kw1/2))*dt;
       Ku2 = a(1)*dt;
       Kw2 = a(2)*dt;  

       [~, xCoor] = min(abs(x_grid(1,:)-(X(i)+(Kx2/2))));
       [~, yCoor] = min(abs(y_grid(:,1)-(Y(i)+(Ky2/2))));
       a = calc_RK_order(U(i), W(i),Ex(yCoor, xCoor),Ey(yCoor, xCoor), q, m);
       
       Kx3 = (U(i)+(Ku2/2))*dt;
       Ky3 = (W(i)+(Kw2/2))*dt;
       Ku3 = a(1)*dt;
       Kw3 = a(2)*dt;

       [~, xCoor] = min(abs(x_grid(1,:)-(X(i)+(Kx3))));
       [~, yCoor] = min(abs(y_grid(:,1)-(Y(i)+(Ky3))));
       a = calc_RK_order(U(i), W(i),Ex(yCoor, xCoor),Ey(yCoor, xCoor), q, m);
       
       Kx4 = (U(i)+Ku3)*dt;
       Ky4 = (W(i)+Kw3)*dt;
       Ku4 = a(1)*dt;
       Kw4 = a(2)*dt; 

       U(i+1) = U(i)+(1/6)*(Ku1+2*Ku2+2*Ku3+Ku4);
       W(i+1) = W(i)+(1/6)*(Kw1+2*Kw2+2*Kw3+Kw4);
       X(i+1) = X(i)+(1/6)*(Kx1+2*Kx2+2*Kx3+Kx4);
       Y(i+1) = Y(i)+(1/6)*(Ky1+2*Ky2+2*Ky3+Ky4);
       
       %----------------------------------- 
       %check if particle hit by electrode.
       %-----------------------------------
       hitElectrode = checkElectrodeProximity( deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ, [Y(i) Y(i+1)], [X(i) X(i+1)], electrodeProximityThresh);
       if (hitElectrode)
           break;
       end
       %-----------------------------------
       
       i=i+1;
   end
end

