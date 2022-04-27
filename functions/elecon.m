function intmel = elecon(nside,intmat,iside,intmel)
% this function fills in the element connectivity matrix needed for the
% fast searching algorithm

    % loop over the sides
    for is = 1:nside % execution loop label 1000
        ip1 = iside(1,is);
        ie1 = iside(3,is);
        ie2 = iside(4,is);
        %first ie1
        i1 = intmat(1,ie1);
        i2 = intmat(2,ie1);
        i3 = intmat(3,ie1);
        if ip1 == i1 
            ipos = 3;
        end
        if ip1 == i2
            ipos = 1;
        end
        if ip1 == i3
            ipos = 2;
        end
        
        % go into intmel
        intmel(ipos,ie1) = ie2;
        % now ie2 if ie2 ~= 0
        if ie2 == 0 % jump gate label 1000
            continue
        else
            i1 = intmat(1,ie2);
            i2 = intmat(2,ie2);
            i3 = intmat(3,ie2); 

            if ip1 == i1 
                ipos = 2;
            end
            if ip1 == i2
                ipos = 3;
            end
            if ip1 == i3
                ipos = 1;
            end
        end
        % go into intmel
        
        intmel(ipos,ie2) = ie1;
        % end loop over the sides
    % jump gate exit label 1000
    end 
    
    return
end