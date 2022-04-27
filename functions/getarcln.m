function [s,ax,ay] = getarcln(t,x0,xc,x1,s)

    % THIS SUBROUTINE COMPUTES THE ARC LENGTH COORDINATE S
    % FOR A BEZIER QUADRATIC SPLINE GIVEN THE VALUE OF THE   
    % T PARAMETER
    temp_x = [];
    % GET AX,BX,AY AND BY
    x = zeros(2,1);

    if ~exist('k','var')
        % third parameter does not exist, so default it to something
        k = 100;
    end

    ax = xc(1) - x0(1);
    bx = x0(1) - 2*xc(1) + x1(1);
    ay = xc(2) - x0(2);
    by = x0(2) - 2*xc(2) + x1(2); 

    a  = bx*bx + by*by;
    b  = 2*(ax*bx + ay*by);
    c  = ax*ax + ay*ay;
    d  = b*b - 4*a*c;



    % IF D=0 THEN IT IS A STRAIGHT LINE 
    if(abs(d) > ((1e-6)*c*c))

        % goto 1000 
        
    else
        % if (k <= 2)
        %     x
        % end
        x = getxy(t,x0,xc,x1);
        temp_x = t;
        s = sqrt((x0(1) - x(1))^2 + (x0(2) - x(2))^2);
        
        return
    end

    % come to 1000

    % GET THE INTEGRATION CONSTANT
    const = (b*sqrt(c))/(4*a);
    const = const - d*log(abs(2*sqrt(a*c) + b))/(8*a*sqrt(a));
    sroot = a*t*t + b*t + c;
    aux   = 2*a*t + b;
    s1    = aux*sqrt(sroot)/(4*a);
    s2    = abs(2*sqrt(a*sroot) + aux);
    s2    = log(s2)/sqrt(a);
    s2    = s2*d/(8*a);

    % NOW GET S BY ADDING UP THE TWO CONTRIBUTIONS

    s = 2*(s1 - s2 - const);

    return
end