function [ G ] = spatial_green( r, z, R_q, Z_q )
%SPATIAL_GREEN Summary of this function goes here
%   Detailed explanation goes here

k =8.987e9;
m = (4*r.*R_q)./(((r+R_q).^2)+((z-Z_q).^2));
G = (k*2./(pi*sqrt((r+R_q).^2 + (z-Z_q).^2)));

r=[]; z =[]; R_q=[]; Z_q =[]; %We clear these bug matrices to save memory

% EDIT NATH: Replaced the commented out lines by the ellipke matlab function:
% simpler to read and more efficient
% teta_num_of_points = 1000;
% teta = linspace(0, pi/2, teta_num_of_points);
% K= zeros(size(G));
% 
% parfor i =1:teta_num_of_points-1
%     diff = 1./(sqrt(1-m*(sin(teta(i)))^2));
%     K = K + diff*pi/(2*(teta_num_of_points-1));
% end
K=ellipke(m);

G = G.* K;

end
