function [ trajParams ] = creatTrajParamsString(q, m, numOfParticles, hitMask, focused, exitR, proxTh, exitRthresh, Ek, entryPr, exitPr, entryR, ...
                                                entryEmittance, exitEmittance)
%CREATTRAJPARAMSSTRING Summary of this function goes here
%   Detailed explanation goes here
c0 = 3e8;
e0 =  -1.60217662e-19;
eM = 9.10938356e-31;
dimFactor = 1e6;
scaleStr = '[\mum]';

exitRRange        = max(abs(exitR(~hitMask)));
entryRRange       = max(abs(entryR(~hitMask)));
particleParams = {' ';
   '-------Particle Parameters-------';
  ['q: ', num2str(q/e0),'[e_0]'];
  ['M: ', num2str(m/eM),'[e_M]'];
  ['Number of Particles: ', num2str(numOfParticles)];
  ' ';};
entryParams ={'----Particle Entry Parameters----';
  ['R_i_n: [', num2str(0),', ',num2str(entryRRange*dimFactor),']',scaleStr];
  ['Ek_i_n: ', num2str(Ek)/1e3,'[KeV]'];
  ['P_r_,_i(max): ',num2str(entryPr)];
  ['\epsilon_i: ', num2str(entryEmittance)];
  ' ';};
exitParams = {'----Particles Exit Parameters----';
  ['R_o_u_t: [', num2str(0),', ',num2str(exitRRange*dimFactor),']',scaleStr];
  ['P_r_,_f(max): ',num2str(exitPr)];
  ['\epsilon_f: ', num2str(exitEmittance)];
  ['d_p_r_o_x: ', num2str(proxTh*dimFactor), scaleStr]
  ['R_f_o_c_u_s: ', num2str(exitRthresh*dimFactor), scaleStr];
  ['Hit: ', num2str(sum(hitMask))];
  ['Focused: ', num2str(focused)];
  };

trajParams = [particleParams; entryParams; exitParams];
end

