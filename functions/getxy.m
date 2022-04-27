function x = getxy(t,x0,xc,x1)

    a1 = (1-t)^2;
    a2 = 2.*t*(1-t);
    a3 = t^2;

    % loop 1000
    for i = 1:2
        x(i) = a1*x0(i)+a2*xc(i)+a3*x1(i);
    end

    return
end