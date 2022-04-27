function xc = intersec(x0,x1,t0,t1,xc,iline)
    
    denom = (t1(2)*t0(1)) - (t1(1)*t0(2));
    
    if (iline == 0)

        % goto 1000
        ;

    else

        % loop 2000
        for j = 1:2
            xc(j) = 0.5*(x0(j)+x1(j)); 
        end
        
        return
    end
    % 1000
    s = ((x1(1)-x0(1))*t0(2)-(x1(2)-x0(2))*t0(1))/denom;

    %3000
    for j = 1:2
        xc(j)=x1(j)+s*t1(j);
    end
    
    return
end