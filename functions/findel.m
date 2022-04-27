function [ie, arZ, i1, i2, i3,ilast] = findel(npoig,neleg,coorg,ieleg,intmeg,x,y,ilast, arZ,i1,i2,i3,ie)
    
    %THIS  FUNCTION   FINDS   OUT  IN  WHICH ELEMENT DOES THE POINT (X,Y) LIE, AND CALCULATES THE AREA COORDINATES FOR INTERPOLATION
    % determinant function
    % start searching with ilast
    imemo = 0;
    if ilast ~= 0 % jump gate label 9

    else
        msg = 'ilast = 0, the background grid needs to be expanded';
        fprintf(msg);
        return % subject to change 
    end

    % jump gate exit label 9
    while true % execution loop label 1000
        ie = ilast;
        i1 = ieleg(1,ie);
        i2 = ieleg(2,ie);
        i3 = ieleg(3,ie);
        x1 = coorg(1,i1);
        x2 = coorg(1,i2);
        x3 = coorg(1,i3);
        y1 = coorg(2,i1);
        y2 = coorg(2,i2);
        y3 = coorg(2,i3);
        %area coordinates
        area2 = deter(x1,y1,x2,y2,x3,y3);
        arZ(1) = deter(x ,y ,x2,y2,x3,y3)/area2;
        arZ(2) = deter(x1,y1,x ,y ,x3,y3)/area2;
        arZ(3) = deter(x1,y1,x2,y2,x ,y )/area2;

        % order them
        for i=1:3 % execution loop label 101
            or(i) = arZ(i);
            ilo(i) = i;
        end

        for i = 1:3 % execution loop label 102
            j1 = i + 1;
            for j = j1:3 % exeuciotn loop label 102, subjected to change
                if (or(i) - or(j) <= 0) % jump gate label 102
                else
                    temp   = or(i);
                    itemp  = ilo(i);
                    or(i)  = or(j);
                    ilo(i) = ilo(j);
                    or(j)  = temp;
                    ilo(j) = itemp;
                end
                % jump gate exit label 102
            end
        end
        
        if (or(1) >= -1*10^-5) % jump gate label 2000
            break
        else
            % get the next element
            inext = intmeg(ilo(1),ie);
            if (inext == 0)
                inext = intmeg(ilo(2),ie);
                if (or(2) < -1*10^-5 && inext ~= 0)
                    ilast = inext;
                    % jump gate label 1000
                    continue
                else

                    if (imemo == 0 ) % jump gate label 880
                    else
                        disp('error in findel')
                        return
                    end     
                end
                % jump gate label 880
                % look for the closest nodal point
                adis = 1*10^6;
                for ip = 1:npoig % execution loop 1010
                    xp = coorg(1,ip);
                    yp = coorg(2,ip);
                    adi1 = (x - xp)^2 + (y - yp)^2;
                    if adi1 > adis % jump gate label 1010
                    else
                        adis = adi1;
                        kpo = ip;
                    end
                    % jump gate exit label 1010
                end  
                % now get an element containing this point
                for ie = 1:neleg % execution loop 1020
                    ilast = ie;
                    
                    for in = 1:3 % execution loop 1021
                
                        ipo = ieleg(in,ie);
                        if(ipo == kpo)
                            pass10 = true; 
                            break
                        else
                            pass10 = false;
                        end
                    end
                    
                    if pass10 == true
                        break
                    end
                end
            else
                ilast = inext;
                continue
            end
        end
        break
    end
    
    % jump gate exit label 2000
    return
end
    