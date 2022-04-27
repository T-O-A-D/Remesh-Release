function [ts,x] = findtsn(nn,x,c,ni,ti,ts)

    % THIS SUBROUTINE COMPUTES THE ARC LENGTH COORDINATE OF THE
    % NTERSECTION POINTS BETWEEN THE BOUNDARY SEGMENT AND 
    % THE BACGROUND GRID

    ii     = single(1);
    ts(ii) = single(0.0);
    sref   = single(0.0);
    s = single([]);

    % loop 1000
    for i = 1:nn-1

        tref  = fix(i-1);
        tref1 = tref+1.;

        % loop 1001
        for j = 1:2 
            x0(j) = x(j,i);
            xc(j) = c(j,i);
            x1(j) = x(j,i+1);
        end

        % loop 1020
        while true
            if ii == ni
                ti(ii+1) = 0;
            end
            if(ti(ii+1) > tref1)

                % goto 1030
                break

            else

                t = ti(ii+1) - tref;
                s = getarcln(t,x0,xc,x1,s);
                ts(ii+1) = sref+s;
                ii = ii+1;
            
                if (ii <= ni)
            
                    % goto 1020
                    continue
                else
                    break

                end

            end

        end

        s    = getarcln(1.0,x0,xc,x1,s);
        sref = sref+s;
    end

    ts(ni) = sref;

    return
end