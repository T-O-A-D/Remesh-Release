
function [npfront,nqfront,iel,nregi] = convex(nelem,npoin,intmat,coord,npfront,nqfront,lcoor,nregi)

    % ADD TRIANGLES TO BACKGROUND GRID TO CREATE CONVEX REGION 

    % FIRST FIND OUT THE BOUNDARY POINTS
    lcoor = boundar(npoin,nelem,lcoor,intmat);
    iyon = 0;
    % STORE BOUNDARY POINTS
    nside = 0;
    for ip = 1:npoin
        if(lcoor(ip) == 0)
            continue
        end
        nside = nside + 1;
        nregi(1,nside) = ip;
    end

    % SET UP THE FRONT
    
    % IDENTIFY A SIDE
    iside = 0;
    iel   = 0;
    while true % execution loop 1999
        iel   = iel + 1;
        ihow  = 0;
        for in = 1:3 % execution loop label 2001
            ip = intmat(in,iel);
            if (lcoor(ip) ~= 0)
                ihow = ihow + 1;
            end
        end
    
        if (ihow >= 2)

        else
            continue
        end

        % GET THE FIRST SIDE      
        for in = 1:3 % execution loop label 3000
            ip = intmat(in,iel);
            if (lcoor(ip) == 0) 
                % jump gate label 3000
            else
                in1 = in + 1;
                if (in1 == 4)
                    in1 = 1;
                end
                ip1 = intmat(in1,iel);
                if (lcoor(ip1) == 0)
                    continue
                    % jump gate label 3000
                else
                    % CHECK FOR THE NEXT FOUR SIDES
                    ip2=ip1+lcoor(ip);
                    if (ip2 <= 0)
                        % jump gate label 3000
                        continue
                    else
                        if (lcoor(ip2) == 0)
                            % jump gate label 3000
                            continue
                        else
                            ip3 = ip + lcoor(ip2);
                            if (ip3 <= 0)
                                % jump gate label 3000
                                continue
                            else
                                if (lcoor(ip3) == 0)
                                    % jump gate label 3000
                                    continue
                                else
                                    ip4 = ip2 + lcoor(ip3);
                                    if (ip4 <= 0)
                                        % jump gate label 3000
                                        continue
                                    else
                                        if (lcoor(ip4) == 0)
                                            % jump gate label 3000
                                            continue
                                        else
                                            ip5 = ip3 + lcoor(ip4);
                                            if (ip5 <= 0)
                                                % jump gate label 3000
                                            else
                                                if (lcoor(ip5) == 0)
                                                    % jump gate label 3000
                                                else
                                                    jump1999 = 0;
                                                    break
                                                    % jump gate label 3001
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            % jump gate exit 3000
            jump1999 = 1;
        end

        if jump1999 == 1
            continue
        end
        % jump gate exit 3001
        iold = ip1;

        % THE FIRST SIDE IS IP1-IP
        % START FILLING THE FRONT
        
        while true % execution loop 4000
            iside = iside + 1;
            npfront(1,iside) = ip1;
            nqfront(1,iside) = ip;
            ikeep = ip;

            if (lcoor(ip) == 0) % jump gate label 4010
                % jump gate exit label 4010
                disp(' ERROR IN CONVEX')
                return
            end
            ip=ip1+lcoor(ip);

            if (ip <= 0) % jump gate label 4010
                % jump gate exit label 4010
                disp(' ERROR IN CONVEX')
                return
            end

            lcoor(ikeep) = 0;
            ip1 = ikeep;

            if (ip1 ~= iold)
                % jump gate label 4000
            else              
                break
            end
        end
        % CHECK IF ISIDE.EQ.NSIDE 
        if (iside == nside) % jump gate label 4500           
            break
        end
    end
    % jump gate exit label 4500
    % LOOP OVER THE NUMBER OF SIDES
    npbou = nside;
    iposi = 0;
    while true % execution loop label 6000
        iposi = iposi+1;
        kn1   = npfront(1,iposi);
        if (kn1 == 0) % jump gate label 5900
            % jump gate exit label 5900
            if (iposi < nside) % jump gate label 6000
                %disp('Jump success')
                continue
            else
                % jump gate label 6000
                return
            end
        end

        kn   = nqfront(1,iposi);
        xn   = coord(1,kn);
        yn   = coord(2,kn);
        xn1  = coord(1,kn1);
        yn1  = coord(2,kn1);
        a    = yn1 - yn;
        b    = xn - xn1;
        c    = (xn1 - xn)*yn1 + (yn - yn1)*xn1;
        inum = 0;
        ang1 = 0.0;
        for is = 1:npbou % execution loop label 6500
            ken = nregi(1,is);
            if (ken == kn || ken == kn1) % jump gate label 6500
                continue
            end
            xken = coord(1,ken);
            yken = coord(2,ken);
            if (a*xken+b*yken+c <= 0.0) % jump gate label 6500
                continue
            end
            xdif1 = xken - xn1;
            ydif1 = yken - yn1;
            xdif2 = xken - xn; 
            ydif2 = yken - yn;
            dist1 = sqrt(xdif1*xdif1+ydif1*ydif1);
            dist2 = sqrt(xdif2*xdif2+ydif2*ydif2);
            cosa  = (xdif1*xdif2 + ydif1*ydif2)/(dist1*dist2);

            if (cosa > 1.0)
                cosa = 1.0;
            end

            if (cosa < -1.0)
                cosa = -1.0;
            end

            angl = acos(cosa);
            if (angl < 0.1) % jump gate label 6500
                continue
            end

            % CHECK WHETHER IT IS A POSSIBLE CONECTIVITY
           
            iyon = possible(kn1,kn,ken,xn1,yn1,xn,yn,xken,yken,1,nside,npbou,npfront,nqfront,nregi,coord,iyon);
            if(iyon == 0)
                % jump gate label 6500
                continue
            end
            inum = 1;
            if (angl < ang1) % jump gate label 6500
                continue
            end
            ang1 = angl;
            kpo = ken;
        end
        if (inum == 0) % jump gate label 5900
            if (iposi < nside)
                % jump gate label 6000
                continue
            else
                return
            end
        end
        % CREATE A NEW ELEMENT
        nelem = nelem + 1;
        intmat(1,nelem) = kn1;
        intmat(2,nelem) = kn;
        intmat(3,nelem) = kpo;

        % UPDATE THE FRONT
        npfront(1,iposi) = 0;
        nqfront(1,iposi) = 0;
        ind = 0;

        for is = 1:nside

            if (npfront(1,is) ~= kpo || nqfront(1,is) ~= kn1)
                break
            end

            ind = 1;
            npfront(1,is) = 0;
            nqfront(1,is) = 0;

            break
        end

        if (ind == 1)
            ind = 0;
        else
            nside = nside + 1; 
            npfront(1,nside) = kn1;
            nqfront(1,nside) = kpo;
        end

        for is = 1:nside

            if (npfront(1,is) ~= kn || nqfront(1,is ~= kpo))
                break
            end

            ind = 1;
            npfront(1,is) = 0;
            nqfront(1,is) = 0;
            break
        end

        if (ind == 1)
            if (iposi < nside)
                % jump gate label 6000
                continue
            else
                return
            end
        else
            break
        end

        nside = nside + 1; 
        npfront(1,nside) = kpo;
        nqfront(1,nside) = kn;
    end
    
    return
end