function [  BetaR, BetaZ, Beta, Gamma ] = getBG( U, W, col )
%return the values of Beta (total and for each axis) and Gamms for all
%particles according to velocity matrices U and W and the pos idx in col;
c0 = 3e8;

row = 1:size(W,1); 
idx = sub2ind(size(W), row, col);

BetaR = W(idx)/c0;
BetaZ = U(idx)/c0;
Beta  = sqrt(BetaZ.^2+BetaR.^2);
Gamma = 1./sqrt(1-Beta.^2);

end

