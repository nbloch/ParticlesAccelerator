function G = analytic_ring_green(r_b, z_b, r_q, z_q)
%Green function: should be multiplied by q to get the potential
k =8.987e9;

m = (4*r_b.*r_q)./(((r_b+r_q).^2)+((z_b-z_q).^2));
% teta_step = 0.001;
% teta = 0:teta_step:pi/2;
% teta= reshape(teta, 1 ,1, length(teta));
% teta= repmat(teta,size(r_b,1), size(r_b,2));
% m = repmat(m,1,1,size(teta,3));
% K = teta_step * sum (1/(sqrt(1 - m.*(sin(teta).^2))),3);
% Was changed to ellipke
K=ellipke(m);
G = (k*2./(pi*sqrt((r_b+r_q).^2 + (z_b-z_q).^2))).* K;

end