function [iin] = interior(iin,x1,y1,x2,y2,x3,y3,x,y)

    % determinant function
    % DETER(P1,Q1,P2,Q2,P3,Q3)=P2*Q3-P3*Q2-P1*Q3+P3*Q1+P1*Q2-P2*Q1

    iin = 1;
    area2 = single(deter(x1,y1,x2,y2,x3,y3));
    if area2 < 1*10^-12

        return

    else

        a1 = single(deter(x,y,x2,y2,x3,y3)/area2);
        a2 = single(deter(x1,y1,x,y,x3,y3)/area2);
        a3 = single(deter(x1,y1,x2,y2,x,y)/area2);
        wcomp = amin1([a1,a2,a3]);

        if wcomp < -0.001 
            iin = 0;
        end
        
        return
    end
end