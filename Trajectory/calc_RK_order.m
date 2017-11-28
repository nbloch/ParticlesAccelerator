function [ a ] = calc_RK_order( Vx,Vy,Ex, Ey, e, m)
%CALC_RK_ORDER Summary of this function goes here
%   Detailed explanation goes here
       c0 = 3e8;
       V_abs = sqrt(Vx^2+Vy^2);
       V_norm = [Vx, Vy]/V_abs;
       V_norm_perp = [V_norm(2), -V_norm(1)];
              
       F = e*[Ex, Ey];
       F_par = F*V_norm.';
       F_perp = F*V_norm_perp.';
       
       gamma = 1/sqrt(1-(V_abs/c0)^2);
       
       a_par  = F_par/((gamma^3)*m);
       a_perp = F_perp/(gamma*m);
       
       a = a_par*V_norm + a_perp*V_norm_perp;
       
end

