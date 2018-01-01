function [  BetaX, BetaZ, BetaY, Beta, Gamma ] = getBG( Vx, Vz, Vy,  col )
%return the values of Beta (total and for each axis) and Gamms for all
%particles according to velocity matrices U and W and the pos idx in col;
c0 = 3e8;

row = 1:size(Vz,1); 
idx = sub2ind(size(Vz), row, col);

BetaX = Vz(idx)/c0;
BetaZ = Vx(idx)/c0;
BetaY = Vy(idx)/c0;
Beta  = sqrt(BetaZ.^2+BetaX.^2+BetaY.^2);
Gamma = 1./sqrt(1-Beta.^2);

end

