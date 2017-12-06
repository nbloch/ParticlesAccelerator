function [ FullParams ] = creatFullParamsString( electrodeWidth, leftElectrodeRadius, rightElectrodeRadius,...
                                               deviceRadius, distanceBetweenElectrodes, deviceLength, ...
                                               VaLeft, VaRight, MSE, M, N, simTraj, trajParams )
%CREATPARAMSSTRING Summary of this function goes here
%   Detailed explanation goes here
DeviceParams = { '--------Device Parameters--------';
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

if(simTraj)
    FullParams = [DeviceParams ; trajParams];
end

end

