function [ hitElectrode ] = checkElectrodeProximity( deviceRadius, leftElectrodeRadius, rightElectrodeRadius, electordesZ, electronPosR, electronPosZ, threshold)
        hitElectrode = 0;
        checkRadius = 0;
        horizDist = abs(electordesZ - electronPosZ(2));
        [M, I] = min(horizDist);
        if(((electronPosZ(2) < electordesZ(I)) && (electronPosZ(1) > electordesZ(I))) || ((electronPosZ(2) > electordesZ(I)) && (electronPosZ(1) < electordesZ(I))))
        checkRadius = 1; %check if electron passed the electrode due to RK low res.
        end
        [~, col, v] = find(horizDist < threshold);
        checkRadius = checkRadius + (~isempty(v));
 if (checkRadius)
      if (mod(col,2) == 1)
            if ((abs(electronPosR(2)) <  deviceRadius) && (abs(electronPosR(2)) >  leftElectrodeRadius-threshold))
                    hitElectrode = 1;
            end
       else
            if ((abs(electronPosR(2)) <  deviceRadius) && (abs(electronPosR(2)) >  rightElectrodeRadius-threshold))
                    hitElectrode = 1;
            end
      end
 end
end
