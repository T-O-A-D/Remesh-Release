function [ttan,amo]= gettan(t,x0,xc,x1,ttan,amo)
    % THIS SUBROUTINE COMPUTES THE POINT TANGENT  COMPONENTS
    % FOR A BEZIER QUADRATIC SPLINE GIVEN THE VALUE OF THE T
    % PARAMETER

    % loop 1000
    for i = 1:2
        ttan(i) = (x0(i) - 2.*xc(i) + x1(i))*t + (xc(i) - x0(i));
    end

    % NORMALIZE
    amo = sqrt(ttan(1)*ttan(1) + ttan(2)*ttan(2));

    % loop 2000
    for i = 1:2
        ttan(i) = ttan(i)/amo;
    end

    return
end