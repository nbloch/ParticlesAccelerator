function [ r_q, z_q , r_b, z_b, leftBoundaryNum, rightBoundaryNum, leftChargesNum, rightChargesNum ] = createUnitCell( M, N, distanceBetweenElectrodes, electrodeWidth, ...
                                                  deviceRadius, leftElectrodeRadius, rightElectrodeRadius )
%CREATEUNITCELL Summary of this function goes here
%   Detailed explanation goes here

%-----------------------------------
%proportions between electrodes:
%-----------------------------------
   leftElectrodeLength  = deviceRadius-leftElectrodeRadius;
   rightElectrodeLength = deviceRadius-rightElectrodeRadius;
   electrodesEffectiveWidth = 2*electrodeWidth;
   leftBoundaryLength  = 2*(electrodesEffectiveWidth + leftElectrodeLength);
   rightBoundaryLength = 2*(electrodesEffectiveWidth + rightElectrodeLength);

   leftChargesWeight  = leftElectrodeLength/(leftElectrodeLength+rightElectrodeLength);
   rightChargesWeight = 1- leftChargesWeight;
   leftBoundaryWeight = leftBoundaryLength/(leftBoundaryLength+rightBoundaryLength);
   rightBoundaryWeight = 1- leftBoundaryWeight;

   leftChargesNum   = round(N*leftChargesWeight);
   rightChargesNum  = N-leftChargesNum;
   leftBoundaryNum  = round(M*leftBoundaryWeight);
   rightBoundaryNum = M-leftBoundaryNum;

%----------------------------
%Charges Spatial Distribution
%----------------------------

RqLeft  = linspace(leftElectrodeRadius  + electrodeWidth, deviceRadius-electrodeWidth, leftChargesNum);
RqRight = linspace(rightElectrodeRadius + electrodeWidth, deviceRadius-electrodeWidth, rightChargesNum);
r_q = [RqLeft, RqRight];

z_q = [(-distanceBetweenElectrodes/2)*ones(1,leftChargesNum), (distanceBetweenElectrodes/2)*ones(1,rightChargesNum)];

%----------------------
%Boundries Distribution
%----------------------

%Left Electrode:
    %relation between num of vertical boundaries and horizontal boundaries
    boundariesHorizProp = 2*electrodesEffectiveWidth/leftBoundaryLength;
    horizBoundNum = round(leftBoundaryNum*boundariesHorizProp);
    %make sure the horizBoundNum will always be a 2 multiplication.
    horizBoundNum = horizBoundNum + mod(horizBoundNum,2);
    %Set at least 2 boundaries condition on the tip of the electrods
    if (horizBoundNum < 2 )
        horizBoundNum = 2;
    end

    vertBoundNum = leftBoundaryNum - horizBoundNum ;

    %Generate vector of vertical boundaries without 2 charges at the corner in
    %order to avoid duplicated charges.
    rVertBoundaryLeft  = linspace(leftElectrodeRadius, deviceRadius, ceil(vertBoundNum/2));
    rVertBoundaryRight = linspace(leftElectrodeRadius, deviceRadius, floor(vertBoundNum/2));

    rHorizBoundaryUp   =  deviceRadius*ones(1,horizBoundNum/2);
    rHorizBoundaryDown =  leftElectrodeRadius*ones(1,horizBoundNum/2);

    zVertBoundaryLeft  = ((-distanceBetweenElectrodes/2) - electrodeWidth)*ones(1,(ceil(vertBoundNum/2)));
    zVertBoundaryRight = ((-distanceBetweenElectrodes/2) + electrodeWidth)*ones(1,(floor(vertBoundNum/2)));
    zHorizBoundaryUp   = linspace((-distanceBetweenElectrodes/2) - electrodeWidth,(-distanceBetweenElectrodes/2) + electrodeWidth, (horizBoundNum/2)+2);
    zHorizBoundaryUp(1)   = [];
    zHorizBoundaryUp(end) = [];
    zHorizBoundaryDown = linspace((-distanceBetweenElectrodes/2) + electrodeWidth,(-distanceBetweenElectrodes/2) - electrodeWidth, (horizBoundNum/2)+2);
    zHorizBoundaryDown(1)   = [];
    zHorizBoundaryDown(end) = [];

    rLeftBoundaryVec = [rVertBoundaryLeft, rHorizBoundaryUp, flip(rVertBoundaryRight), rHorizBoundaryDown];
    zLeftBoundaryVec = [zVertBoundaryLeft, zHorizBoundaryUp, zVertBoundaryRight, zHorizBoundaryDown];

%Right Electrode:
    %relation between num of vertical boundaries and horizontal boundaries
    boundariesHorizProp = 2*electrodesEffectiveWidth/rightBoundaryLength;
    horizBoundNum = round(rightBoundaryNum*boundariesHorizProp);
    %make sure the horizBoundNum will always be a 2 multiplication.
    horizBoundNum = horizBoundNum + mod(horizBoundNum,2);
    %Set at least 2 boundaries condition on the tip of the electrods
    if (horizBoundNum < 2 )
        horizBoundNum = 2;
    end

    vertBoundNum = rightBoundaryNum - horizBoundNum ;

    %Generate vector of vertical boundaries without 2 charges at the corner in
    %order to avoid duplicated charges.
    rVertBoundaryLeft = linspace(rightElectrodeRadius, deviceRadius, ceil(vertBoundNum/2));
    rVertBoundaryRight = linspace(rightElectrodeRadius, deviceRadius, floor(vertBoundNum/2));

    rHorizBoundaryUp   =  deviceRadius*ones(1,horizBoundNum/2);
    rHorizBoundaryDown =  rightElectrodeRadius*ones(1,horizBoundNum/2);

    zVertBoundaryLeft  = ((distanceBetweenElectrodes/2) - electrodeWidth)*ones(1,(ceil(vertBoundNum/2)));
    zVertBoundaryRight = ((distanceBetweenElectrodes/2) + electrodeWidth)*ones(1,(floor(vertBoundNum/2)));
    zHorizBoundaryUp   = linspace((distanceBetweenElectrodes/2) - electrodeWidth,(distanceBetweenElectrodes/2) + electrodeWidth, (horizBoundNum/2)+2);
    zHorizBoundaryUp(1)   = [];
    zHorizBoundaryUp(end) = [];
    zHorizBoundaryDown = linspace((distanceBetweenElectrodes/2) + electrodeWidth,(distanceBetweenElectrodes/2) - electrodeWidth, (horizBoundNum/2)+2);
    zHorizBoundaryDown(1)   = [];
    zHorizBoundaryDown(end) = [];

    rRightBoundaryVec = [rVertBoundaryLeft, rHorizBoundaryUp, flip(rVertBoundaryRight), rHorizBoundaryDown];
    zRightBoundaryVec = [zVertBoundaryLeft, zHorizBoundaryUp, zVertBoundaryRight, zHorizBoundaryDown];


r_b = [rLeftBoundaryVec, rRightBoundaryVec].';
z_b = [zLeftBoundaryVec, zRightBoundaryVec].';

%----------------------
%Summary
%----------------------
r_q = repmat(r_q, M, 1);
z_q = repmat(z_q, M, 1);
r_b = repmat(r_b, 1, N);
z_b = repmat(z_b, 1, N);

end
