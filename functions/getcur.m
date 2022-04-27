function [curv,amo,ttan] = getcur(t,x0,xc,x1,curv,amo,ttan)

    % THIS  SUBROUTINE  COMPUTES   THE   POINT   CURVATURE
    % FOR A BEZIER QUADRATIC SPLINE GIVEN THE VALUE OF THE
    % T PARAMETER

    % GET THE TANGENT
    [ttan,amo] = gettan(t,x0,xc,x1,ttan,amo);

    % THE NORMAL
    anr(1) =  ttan(2);
    anr(2) = -ttan(1);
    curv   = 0.0;

    % loop 1000
    for i = 1:2
        curv = curv + (x0(i) - 2.*xc(i) + x1(i))*anr(i);
    end

    curv = abs(2.*curv*amo*amo);

    return
end