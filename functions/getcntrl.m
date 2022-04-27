
function [x,c] = getcntrl(n,x,c,lx,t,ls, xc)

    
    % THIS SUBROUTINE EVALUATES THE CONTROL POINTS TO OBTAIN A
    % CUADRATIC BEZIERS FITTING.
    ls = rfilliv(ls,n,0);

    % 
    if (n == 2)
        % goto 6800
        % GET FINAL CONTROL POINTS
        for i =1:n-1
            for j =1:2
                x0(j) = x(j,i);
                x1(j) = x(j,i+1);
                t0(j) = t(j,i);
                t1(j) = t(j,i+1);
            end

            % FIND INTERSECTION
            iline = 0;

            if(ls(i) == 0 && ls(i+1) == 0)
                iline = 1;
            end
        
            xc = intersec(x0,x1,t0,t1,xc,iline);
            % loop 6002
            for j =1:2
                c(j,i) = xc(j);
            end   
        end
        
        return
    end

    % IDENTIFY STRAIGHT LINE SEGMENTS
    % loop 500
    for i = 2:n-1
        % loop 501
        for j = 1:2
            t0(j) = x(j,i)   - x(j,i-1);
            t1(j) = x(j,i+1) - x(j,i);
        end
        % SCALAR PRODUCT
        sca = t0(1)*t1(1) + t0(2)*t1(2);

        if (sca <= 0.0)
            % loop 5000
            disp("' *** MORE POINTS ARE NEEDED TO DEFINE THE BOUNDARY'")
            return
        end

        % FIND MODULUS
        am0  = sqrt(t0(1)*t0(1) + t0(2)*t0(2)); 
        am1  = sqrt(t1(1)*t1(1) + t1(2)*t1(2)); 
        vect = (t0(2)*t1(1)-t0(1)*t1(2))/(am0*am1);

        if(abs(vect) < 1.e-3) 
            %goto 500
            continue
        end

        if(vect > 1.e-3)
            ls(i) = 1;
        end

        if(vect < -1.e-3) 
            ls(i) = -1;
        end
    end

    ls(1) = 2;
    ls(n) = 2;

    lx = rfilliv(lx,n,0);

    % LOOP OVER SEGMENTS
    if (n == 3)
        % goto 601
        ;
    else
        % loop 600
        for i = 2:n-2

            % IF THERE IS CHANGE IN CURVATURE, MARK THEM
            ilo = ls(i)*ls(i+1);

            if (ilo < 0)
                lx(i) = 1;
            end

        end

        % CHECK NEIGHBOURS
        % loop 700
        for i = 2:n-2

            if(lx(i) == 0) 
                continue
            end
            
            if(lx(i-1) == 1 || lx(i+1) == 1) 
                disp("' *** MORE POINTS ARE NEEDED TO DEFINE THE BOUNDARY'")
                return
            end

            if(ls(i-1) == 0 || ls(i+1) == 0) 
                disp("' *** MORE POINTS ARE NEEDED TO DEFINE THE BOUNDARY'")
                return
            end
        end
    end

    % THE POINTS ARE ALL RIGHT !!!!
    % GET FIRST THE TANGENTS PARALEL TO THE SIDES 
    % loop 680
    for i = 1:n-1

        if (ls(i) == 0)
            continue
        end

        if (ls(i+1) ~= 0)
            continue
        end

        ls(i) = 0;

        % loop 690
        for j = 1:2
            t(j,i) = x(j,i+1) - x(j,i);
        end
    end

    % loop 681
    for k = 1:n-1

        i = n + 1 - k;

        if (ls(i) == 0)
            continue
        end

        if (ls(i+1) ~= 0)
            continue
        end

        ls(i) = 0;

        % loop 691
        for j = 1:2
            t(j,i) = x(j,i) - x(j,i-1);
        end
    end

    % loop 800
    for i = 2:n-2

        if (lx(i) == 0)
            continue
        end

        % loop 801
        for j = 1:2
            a(j) = x(j,i+1) - x(j,i);
        end

        % loop 802
        for k = i:i+1
            ls(k) = 0;

            % loop 803
            for j = 1:2
                t(j,k) = a(j);
            end
        end
    end

    % FIND OUT THE POINT SLOPES

    % loop 1000
    for i = 1:n-1

        % loop 1001
        for j = 1:2
            a(j) = x(j,i+1) - x(j,i);
        end

        % ADD THEM

        % loop 1002
        for k = i:i+1
            if (ls(k) == 0)
                continue
            end

            % loop 1003
            for j = 1:2
                t(j,k) = t(j,k) + a(j);
            end
        end
    end

    % GET THE END SLOPES
    % NO LINK BETWEEN 1 AND N ----> GET TANGENT VECTORS

    % loop 3101
    for j = 1:2
        a(j) = x(j,2) - x(j,1);
    end

    amodu = sqrt(a(1)*a(1) + a(2)*a(2));

    % loop 3102
    for j = 1:2
        a(j) = a(j)/amodu;
    end

    sca = 0.0;

    % loop 3103
    for j = 1:2
        sca = sca + a(j)*t(j,2);
    end

    % loop 3104
    for j = 1:2
        t(j,1) = t(j,2) + 2.*(sca*a(j) - t(j,2));
    end

    % loop 3111
    for j = 1:2
        a(j) = x(j,n) - x(j,n-1);
    end

    amodu = sqrt(a(1)*a(1) + a(2)*a(2));

    % loop 3112
    for j = 1:2
        a(j) = a(j)/amodu;
    end

    sca = 0.0;

    % loop 3113
    for j = 1:2
        sca = sca + a(j)*t(j,n-1);
    end

    % loop 3114
    for j = 1:2
        t(j,n) = t(j,n-1) + 2.*(sca*a(j) - t(j,n-1));
    end
    
    % GET FINAL CONTROL POINTS
    % loop 6000
    for i =1:n-1
        % loop 6001
        for j =1:2
            x0(j) = x(j,i);
            x1(j) = x(j,i+1);
            t0(j) = t(j,i);
            t1(j) = t(j,i+1);
        end

        % FIND INTERSECTION
        iline = 0;

        if(ls(i) == 0 && ls(i+1) == 0)
            iline = 1;
        end

        xc = intersec(x0,x1,t0,t1,xc,iline);

        % loop 6002
        for j =1:2
            c(j,i) = xc(j);
        end

    end

    return
end