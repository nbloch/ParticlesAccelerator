function [ trajParams ] = creatTrajParamsString(q, m, numOfParticles, hitMask, focused, endBetaZ,...
                                                endBetaR, exitR, proxTh, exitRthresh, lowAxialVel,...
                                                highAxialVel, lowRadialVel,highRadialVel, entryR)
%CREATTRAJPARAMSSTRING Summary of this function goes here
%   Detailed explanation goes here
c0 = 3e8;
e0 =  -1.60217662e-19;
eM = 9.10938356e-31;

exitLowAxialVel   = min(endBetaZ(~hitMask));
exitHighAxialVel  = max(endBetaZ(~hitMask));
exitLowRadialVel  = min(abs(endBetaR(~hitMask)));
exitHighRadialVel = max(abs(endBetaR(~hitMask)));
exitRRange        = max(abs(exitR(~hitMask)));
entryRRange       = max(abs(exitR(~entryR)));
trajParams = {' ';
   '-------Particle Parameters-------';
  ['q: ', num2str(q/e0),'[e_0]'];
  ['M: ', num2str(m/eM),'[e_M]'];
  ['Number of Particles: ', num2str(numOfParticles)];
  ' ';
   '----Particle Entry Parameters----';
  ['V_z_-_i_n: [', num2str(lowAxialVel/c0),', ',num2str(highAxialVel/c0),'][c]'];
  ['V_r_-_i_n: [', num2str(lowRadialVel/c0),', ',num2str(highRadialVel/c0),'][c]'];
  ['R_i_n: [', num2str(-entryRRange*1e3),', ',num2str(entryRRange*1e3),'][mm]'];
  ' ';
   '----Particles Exit Parameters----';
  ['V_z_-_o_u_t: [', num2str(exitLowAxialVel),', ',num2str(exitHighAxialVel),'][c]'];
  ['V_r_-_o_u_t: [', num2str(exitLowRadialVel),', ',num2str(exitHighRadialVel),'][c]'];
  ['R_o_u_t: [', num2str(-exitRRange*1e3),', ',num2str(exitRRange*1e3),'][mm]'];
  ['d_p_r_o_x: ', num2str(proxTh*1e6), '[\mum]']
  ['R_f_o_c_u_s: ', num2str(exitRthresh*1e3), '[mm]'];
  ['Hit: ', num2str(sum(hitMask))];
  ['Focused: ', num2str(focused)];
  };

end

