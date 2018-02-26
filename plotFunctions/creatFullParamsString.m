function [ DeviceParams ] = creatFullParamsString( electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                               deviceRadius, distanceBetweenElectrodes, deviceLength, ...
                                               VaLeft, VaRight, MSE, M, N)
%CREATPARAMSSTRING Summary of this function goes here
%   Detailed explanation goes here
scaleFactor = 1e6;
scaleStr ='[\mum]' ;
DeviceParams = { '--------Device Parameters--------';
              ['Electrode Width: ', num2str(electrodeWidth*scaleFactor),scaleStr];
              ['R_a_p_-_l: ', num2str(leftElectrodeRadius*scaleFactor),scaleStr];
              ['R_a_p_-_r: ', num2str(rightElectrodeRadius*scaleFactor),scaleStr];
              ['R_L_e_n_s: ', num2str(deviceRadius*scaleFactor),scaleStr];
              ['d_a_p: ', num2str(distanceBetweenElectrodes*scaleFactor),scaleStr];
              ['l_L_e_n_s: ', num2str(deviceLength*scaleFactor),scaleStr];
               '-------Electric Parameters-------';
              ['V_a_-_L_e_f_t: ', num2str(VaLeft),'[V]'];
              ['V_a_-_R_i_g_h_t: ', num2str(VaRight),'[V]'];
               '------Resolution Parameters------';
              ['MSE: ', num2str(MSE),'[%]'];
              ['M_b: ', num2str(M)];
              ['N_c: ', num2str(N)];
		};

end

