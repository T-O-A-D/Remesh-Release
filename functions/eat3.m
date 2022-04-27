function [lcore,icone,coord,iside,lcoid,intmat,nelem,iloca] = eat3(npoin,nelem,nside,intmat,iside,lcoor,lposi,lwher,lhowm,lcore,lcoid,icone,coord)
    % this function removes the points where only three elements coincide
    ntres = 0;
    kpoin = 0;

    % loop over number of nodes
    for ip = 1:npoin % execution loop label 1000

        kpoin = kpoin + 1;
        lposi(ip) = kpoin;

        % if boundary point leave it
        if lcoor(ip) ~= 0 % jump gate label 1000
            continue
        else
            % check number of elements
            if lcore(ip) ~= 3 % jump gate label 1000
                continue
            else
                % a point with only three elements
                % get the elements from icone
                iloca = lwher(ip);
                ie1   = icone(iloca+1);
                ie2   = icone(iloca+2);
                ie3   = icone(iloca+3);
                if (ie1 == 0 || ie2 == 0 || ie3 == 0) % jump gate label 1000
                    continue
                else
                    ntres = ntres + 1;
                    kpoin = kpoin - 1;
                    lposi(ip) = 0;
                    % get new connectivity point from ie1 from ie2
                    ip1 = intmat(1,ie1);
                    ip2 = intmat(2,ie1);
                    ip3 = intmat(3,ie1);
                    for in = 1:3 % execution loop label 1002
                        ipt = intmat(in,ie2);
                        if (ipt ~= ip1 && ipt ~= ip2 && ipt ~= ip3)
                            ino = ipt;
                        end
                    end % execution loop label 1002
                    
                    % replace connectivity
                    for in = 1:3 % execution loop label 1003
                        if (intmat(in,ie1) == ip )
                            intmat(in,ie1) = ino;
                        end
                        lcore(intmat(in,ie1)) = lcore(intmat(in,ie1)) - 1;
                    end
                    
                    for in = 1:3 % execution loop label 1004
                        intmat(in,ie2) = 0;
                        intmat(in,ie3) = 0;
                    end
                    
                    % end loop over points 
                    iloca = lwher(ino);
                    for in = 1:lcore(ino) % execution loop label 1005
                        if ((icone(iloca+in) == ie2) || (icone(iloca+in) == ie3) )
                            icone(iloca+in) = 0;
                        end
                    end
                end       
            end               
        end
    end

    % jump gate exit label 1000
    % transfer
    for ip = 1:npoin % execution loop label 2000
        il = lposi(ip);
        if (il == 0) % jump gate label 2000
            continue
        else
            lcoor(il) = lcoor(ip);
            lcore(il) = lcore(ip);
            lcoid(il) = lcoid(ip);
            for id = 1:2 % execution loop label 2001
                coord(id,il) = coord(id,ip);
            end
        end
        % jump gate exit label 2000
    end

    je = 0;
    for ie = 1:nelem % execution loop label 3000
        if intmat(1,ie) == 0 
            continue
        end
        je = je + 1;
        for in = 1:3 % execution loop label 3001
            iold = intmat(in,ie);
            inew = lposi(iold);
            intmat(in,je) = inew;
        end
    end
    % get npoin and nelem
    npoin = npoin - ntres;
    nelem = nelem - 2*ntres;

    %output nr. of eaten points
    fprintf('nr. of 3''s removed = %d\n',ntres)    
    % fill in iside
    if ntres > 0 
        fprintf('Filling iside\n')
        [iside,iloca,lwher,icone] = side(nelem,npoin,nside,intmat,iside,lwher,lhowm,icone);
    else
        iloca = nside;
    end

    return
end