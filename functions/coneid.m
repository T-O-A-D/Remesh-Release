function lcoid = coneid(npoin,nelem,intmat,lcoor,lcoid,coord,strec)

    % This sub find out the optimal number of connectivities for each mode

    piZ = 3.1415927;
    piZ = single(piZ);

    % loop 1000
    for ip = 1:npoin
        lcoid(ip) = 6;
    end

    % searc for a side to start
    for ip = 1:npoin % execution loop label 2000

        ic = ip;

        if lcoor(ip) ~= 0 % jump gate label 2001
            break
        end
    end

    % jump gate label 2001
    ia = 0;
    ib = 0;
    superbreak = 0;

    for ie = 1:nelem % execution loop label 3000
        for in = 1:3 % execution loop label 3001

            ip = intmat(in,ie);

            if ip ~= ic % jump gate label 3001

                continue;
            end
            
            inb = in - 1;

            if inb < 1 
                inb = inb + 3;
            end

            if lcoor(intmat(inb,ie)) ~= 0
                ib = intmat(inb,ie);
            end

            if ib ~= 0 % jump gate label 3002
                superbreak = 1;
                break
            end

            
            % jump gate exit label 3001
        end

        if superbreak == 1
            break
        end

    end

    % jump gate exit label 3002
    imemo = ic;

    while true

        ia   = ib-lcoor(ic);
        alph = strec(2,ic);
        anx  = strec(3,ic);
        any  = strec(4,ic);
        x1r  = coord(1,ib)-coord(1,ic);
        y1r  = coord(2,ib)-coord(2,ic);
        x2r  = coord(1,ia)-coord(1,ic);
        y2r  = coord(2,ia)-coord(2,ic);
        x1   = fxba(x1r,y1r,anx,any,alph);
        y1   = fyba(x1r,y1r,anx,any,alph);
        x2   = fxba(x2r,y2r,anx,any,alph);
        y2   = fyba(x2r,y2r,anx,any,alph);
        cosa = (x1*x2+y1*y2)/(sqrt(x1*x1+y1*y1)*sqrt(x2*x2+y2*y2));

        if cosa > 1.0
            cosa = 1;
        end

        if cosa < -1.0
            cosa = -1;
        end
        
        theta = acos(cosa);

        if (x2*y1 - x1*y2) < 0
            theta = 2*piZ - theta;
        end

        divi = floor((theta + (piZ/6))/(piZ/3));
        
        lcoid(ic) = divi;

        if lcoid(ic) < 1
            lcoid(ic) = 1;
        end

        ib = ic;
        ic = ia;

        if ic ~= imemo
            continue
        else
            break;
        end
    end
    
    return
end