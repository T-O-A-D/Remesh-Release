
function coorg = expand(neleg,npoig,coorg,coor0,ieleg,intmeg,n1body,n2body,unit_7,fid_7,cv)
    
    coor0 = zeros(2,npoig);
    % MODIFY COORDINATES OF THE BACKGROUND GRID BOUNDARY POINTS TO 
    % COMPLETELY COVER THE DOMAIN OF INTEREST
    % get EXPA
    
    expa = input('Input expansion factor (generally close to 0) :\n');
    
    % -------  to single -------
    expa = single(expa);
    
    coor0 = rfillm(coor0,2,npoig,0.0);
    
    % LOOP OVER THE ELEMENTS
    ae = 1.e+6;

    % loop 1000
    for ie = 1:neleg
        % loop 1001
        for in = 1:3
            if (intmeg(in,ie) ~= 0)
                % goto 1001
                continue
            end

            % MOVE NODES
            in1 = in + 1;
            in2 = in + 2;

            if (in1 > 3)
                in1 = in1 - 3;
            end

            if (in2 > 3)
                in2 = in2 - 3;
            end

            ip1 = ieleg(in1,ie);
            ip2 = ieleg(in2,ie);
            ax  = coorg(2,ip2) - coorg(2,ip1);
            ay  = coorg(1,ip1) - coorg(1,ip2);
            at  = sqrt(ax*ax + ay*ay);
            
            if (at < ae)
                ae =at;
            end
        end
    end
    disp('AE = ')
    disp(ae)
    
    % loop 2000
    for ie = 1:neleg
        exp = expa*ae;

        if (ie >= n1body && ie <= n2body)
            exp = exp*1.e-2;
        end
        % loop 2001
        for in = 1:3

            if (intmeg(in,ie) ~= 0)
                continue
            end

            % MOVE NODES
            in1 = in + 1;
            in2 = in + 2;

            if (in1 > 3)
                in1 = in1 - 3;
            end

            if (in2 > 3)
                in2 = in2 - 3;
            end

            ip1          = ieleg(in1,ie);
            ip2          = ieleg(in2,ie);
            ax           = coorg(2,ip2) - coorg(2,ip1);
            ay           = coorg(1,ip1) - coorg(1,ip2);
            at           = sqrt(ax*ax + ay*ay); 
            ax           = ax*exp/at;
            ay           = ay*exp/at;
            coor0(1,ip1) = coor0(1,ip1) + ax;
            coor0(2,ip1) = coor0(2,ip1) + ay;
            coor0(1,ip2) = coor0(1,ip2) + ax;
            coor0(2,ip2) = coor0(2,ip2) + ay;
        
        end
    end

    % ADDING
    % 3000
    for i =1:npoig
        % 3001
        for j =1:2

            coorg(j,i) = coorg(j,i) + coor0(j,i);

        end
    end

    return
end