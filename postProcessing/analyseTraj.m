function [trajLen, hit, focused, hitMask, focusedMask ] = analyseTraj( X, Y, Zq, deviceRadius, leftElectrodeRadius,...
                                        rightElectrodeRadius, electrodeProximityThresh, exitRthresh, lastZ, numOfParticles)
%ANALYSETRAJ Summary of this function goes here
%   Detailed explanation goes here

trajLen = (sum(~isnan(X),2))';
electordesZ = unique(Zq);
focused = 0;
focusedMask = zeros(1,numOfParticles);

parfor i=1:numOfParticles
    %check If Particle i hit the electrode
    hitMask(i) = checkElectrodeProximity( deviceRadius, leftElectrodeRadius, rightElectrodeRadius,...
                                          electordesZ, [Y(i,trajLen(i)-1) Y(i,trajLen(i))],...
                                          [X(i,trajLen(i)-1) X(i,trajLen(i))],electrodeProximityThresh);
    %check If Particle i passed through the focuse radius
    if ((abs(Y(i,trajLen(i))) < exitRthresh) && (X(i,trajLen(i)) >= lastZ))
        focused = focused +1;
        focusedMask(i) = 1;
    end
end

hit = sum(hitMask);
passMask = ~hitMask;
notFocusedMask = logical(passMask-focusedMask);
focusedMask = logical(focusedMask);

end

