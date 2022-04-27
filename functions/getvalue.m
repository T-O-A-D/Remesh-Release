function [dis, alp, anx, any] = getvalue(npoig,i1,i2,i3,arZ,delta,dis,alp,anx,any,inorm)

    % interpolate
    dis = single(arZ(1)*delta(1,i1)+arZ(2)*delta(1,i2)+arZ(3)*delta(1,i3));
    alp = single(arZ(1)*delta(2,i1)+arZ(2)*delta(2,i2)+arZ(3)*delta(2,i3));
    anx = single(arZ(1)*delta(3,i1)+arZ(2)*delta(3,i2)+arZ(3)*delta(3,i3));
    any = single(arZ(1)*delta(4,i1)+arZ(2)*delta(4,i2)+arZ(3)*delta(4,i3));


    %normalize if required
    if(inorm == 0) 
        return
    end
    anm = sqrt(anx.*anx + any.*any);

    anx = anx/anm;
    any = any/anm;
    
    return
end 
