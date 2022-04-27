function iyon = possible(kn1,kn,kp,xn1,yn1,xn,yn,xp,yp,ireg,nonf,nonr,npfront,nqfront,nregi,coor,iyon)

    % This function finds out whether connection with point kp is possible iyon = 1 or not iyon = 0 %
    iin = [];
    iyon = 1;
    kn1 = single(kn1);
    kn = single(kn);
    kp = single(kp);
    xn1 = single(xn1);
    yn1 = single(yn1);
    xn = single(xn);
    yn = single(yn);
    xp = single(xp);
    yp = single(yp);
    ireg = single(ireg);
    nonf = single(nonf);
    nonr = single(nonr);
    npfront = single(npfront);
    nqfront = single(nqfront);
    nregi = single(nregi);
    coor = single(coor);
    % loop over the front nodes
    for it = 1:nonr(ireg) %execution loop label 1000

        kj = single(nregi(ireg,it));

        if (kj == kn1 || kj == kn || kj == kp) % jump gate label 1000
            continue
        else

            xt = coor(1,kj);
            yt = coor(2,kj);

            % check if the point is interior
            iin = single(interior(iin,xn1,yn1,xn,yn,xp,yp,xt,yt));

            if iin == 0 % jump gate label 1000
                continue
            else    
                iyon = 0;
                return
            end
        end 
    % jump gate exit label 1000 
    end

    % equation of the mid-base : kp line
    xmb = 0.5*(xn1+xn);
    ymb = 0.5*(yn1+yn);
    as  = ymb-yp;
    bs  = xp-xmb;
    cs  = (xmb-xp)*ymb+(yp-ymb)*xmb;

    % loop over the front sides : check for intersection

    for ir = 1:nonf(ireg) % execution loop label 2000
        
        knt1 = npfront(ireg,ir);

        if knt1 == 0 % jump gate label 2000
            continue
        else

            knt = nqfront(ireg,ir);

            if (knt1 == kn1 && knt == kn) % jump gate label 2000
                continue
            else
                
                if (knt1 == kp || knt == kp) % jump gate label 2000
                    continue
                else

                    xnt1 = coor(1,knt1);
                    ynt1 = coor(2,knt1);
                    xnt  = coor(1,knt);
                    ynt  = coor(2,knt);
                    at   = ynt1-ynt;
                    bt   = xnt-xnt1;
                    ct   = (xnt1-xnt)*ynt1+(ynt-ynt1)*xnt1;
                    s1   = at*xmb+bt*ymb+ct;
                    s2   = at*xp+bt*yp+ct;
                    s3   = as*xnt1+bs*ynt1+cs;
                    s4   = as*xnt+bs*ynt+cs;
                    sig1 = s1*s2;
                    sig2 = s3*s4;

                    if (sig1 > 0.0 || sig2 > 0.0) % jump gate label 2000
                        continue
                    else
                        iyon = 0;
                        return
                    end
                end
            end
        end
    % jump gate exit label 2000 
    end
    
    return
end