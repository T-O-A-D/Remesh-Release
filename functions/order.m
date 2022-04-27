function [npfront, nqfront, disw] = order(nl,ireg,npfront,nqfront,coor,npoig,neleg,coorg,ieleg,delta,strec,disw,disw1) 

    disw = 1.e+6;
    coor = single(coor);
    for il = 1:nl % execution loop label 1000
        
        kn    = single(nqfront(ireg,il));
        kn1   = single(npfront(ireg,il));
        
        xn    = coor(1,kn) ;
        yn    = coor(2,kn) ;
        xn1   = coor(1,kn1);
        yn1   = coor(2,kn1);
        
        alph  = single(amin1([(strec(2,kn1)),strec(2,kn)]));
        anx   = single(0.5*(  strec(3,kn1)+strec(3,kn)) );
        any   = single(0.5*(  strec(4,kn1)+strec(4,kn)) );

        anm   = sqrt(anx*anx+any*any);
        anx   = anx/anm;
        any   = any/anm;

        xnf2  = fxba(xn,yn,anx,any,alph);
        ynf2  = fyba(xn,yn,anx,any,alph);
        xn1f2 = fxba(xn1,yn1,anx,any,alph);
        yn1f2 = fyba(xn1,yn1,anx,any,alph);

        x12   = xnf2-xn1f2;
        y12   = ynf2-yn1f2;
        dis   = sqrt(x12^2+y12^2);
        

        if dis < disw1 % jump gate label 1000
            
            continue;
        end

        if dis > disw % jump gate label 1000
            
            continue;
        end

        iwh  = il;
        disw = dis;
        iq   = kn;
        ip   = kn1;
        
    end

    % jump gate exit label 1000
    % swap values

    npfront(ireg,iwh) = npfront(ireg,nl);
    nqfront(ireg,iwh) = nqfront(ireg,nl);
    npfront(ireg,nl)  = ip;
    nqfront(ireg,nl)  = iq;

    return
end