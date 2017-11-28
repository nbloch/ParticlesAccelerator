function    [r_q, z_q, r_b, z_b, z_off] = createLastElectrode(M, N, distanceBetweenElectrodes, electrodeWidth, ...
                                                  deviceRadius, electrodeRadius)
%CREATELASTELECTRODE Summary of this function goes here
%   Detailed explanation goes here

%-----------------------------------
%proportions between electrodes:
%-----------------------------------

electrodeLength  = deviceRadius-electrodeRadius;
electrodesEffectiveWidth = 2*electrodeWidth;
boundaryLength  = 2*(electrodesEffectiveWidth + electrodeLength);

%----------------------------
%Charges Spatial Distribution
%----------------------------

r_q  = linspace(electrodeRadius  + electrodeWidth, deviceRadius-electrodeWidth, N);
z_q = zeros(1,N);

%----------------------
%Boundries Distribution
%----------------------

boundariesHorizProp = 2*electrodesEffectiveWidth/boundaryLength;
horizBoundNum = round(M*boundariesHorizProp);
%make sure the horizBoundNum will always be a 2 multiplication.
horizBoundNum = horizBoundNum + mod(horizBoundNum,2);
%Set at least 2 boundaries condition on the tip of the electrods
if (horizBoundNum < 2 )  
    horizBoundNum = 2;
end

vertBoundNum = M - horizBoundNum ;

%Generate vector of vertical boundaries without 2 charges at the corner in
%order to avoid duplicated charges.
rVertBoundaryLeft  = linspace(electrodeRadius, deviceRadius, ceil(vertBoundNum/2));
rVertBoundaryRight = linspace(electrodeRadius, deviceRadius, floor(vertBoundNum/2));

rHorizBoundaryUp   =  deviceRadius*ones(1,horizBoundNum/2);
rHorizBoundaryDown =  electrodeRadius*ones(1,horizBoundNum/2);

zVertBoundaryLeft  = -electrodeWidth*ones(1,(ceil(vertBoundNum/2)));
zVertBoundaryRight =  electrodeWidth*ones(1,(floor(vertBoundNum/2)));
zHorizBoundaryUp   = linspace(-electrodeWidth,  electrodeWidth, (horizBoundNum/2)+2);
zHorizBoundaryUp(1)   = [];
zHorizBoundaryUp(end) = [];
zHorizBoundaryDown = linspace( electrodeWidth, -electrodeWidth, (horizBoundNum/2)+2);
zHorizBoundaryDown(1)   = [];
zHorizBoundaryDown(end) = [];    

r_b = [rVertBoundaryLeft, rHorizBoundaryUp, flip(rVertBoundaryRight), rHorizBoundaryDown].';  
z_b = [zVertBoundaryLeft, zHorizBoundaryUp, zVertBoundaryRight, zHorizBoundaryDown].';

z_off = 1.5*distanceBetweenElectrodes;
%z_off is relative to the last unit cell center

%----------------------
%Summary
%----------------------
r_q = repmat(r_q, M, 1);
z_q = repmat(z_q, M, 1);
r_b = repmat(r_b, 1, N);    
z_b = repmat(z_b, 1, N);

end

