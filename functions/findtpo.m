function [tx, xx, np, tl, srein,tg, x, ax, ay] = findtpo(nn,ni,np,x,c,ts,tg,td,tx,xx ,tl, srein)

    %this subroutine computes the  number and the t coordinate  of  the
    %     generated boundary points.  the process is  as follows ; integrate
    %     the spacing: tg ---> determine number of points: np ---> determine
    %     their s coordinate: tl ---> determine their t coordinate: tx  --->
    %     determine the x and y coordinates: xx.

    temp_cres = [];
    temp_x = [];
    cres = single(0);
    niter = single(15);
    % integrate
 
    tg(1) = single(0.0);

    for ii = 2:ni % execution loop label 1000

        sint   = ts(ii) - ts(ii-1);
        tg(ii) = tg(ii-1) + 0.5*(td(ii-1) + td(ii))*sint;
        
    end
    tg = tg(1:ni);
    ts = ts(1:ni);
    

    % find out number of points
    np = fix(tg(ni)) + 1;
    sp = tg(ni)/real(np);
    np = np +1;

    % find the s coordinate
    ip = 1;
    tl(ip) = 0.0;
    case2001 = true;

    for ii = 2:ni % execution loop label 2000
        
        while case2001 == true
            
            ts1 = real(ip)*sp;
            
            if (tg(ii) < ts1) % % jump gate label 2000
                
                break

            else

                sl   = ts1 - tg(ii-1);
                a1   = 0.5*(td(ii) - td(ii-1))/(ts(ii) - ts(ii-1));
                sint = 0.5*(ts(ii) - ts(ii-1));
                stry = sint;

                for k = 1:niter % execution loop label 5000

                    sint = 0.5*sint;
                    sres = td(ii-1)*stry+a1*stry^2;

                    if sres-sl <= 0

                        stry = stry + sint;

                    else

                        stry = stry - sint;

                    end
                end

                ip = ip + 1;
                
                tl(ip) = ts(ii-1)+stry;
                % jump gate label 2001

            end
        end
        % jump gate exit label 2000
    end

    
    tl(np) = ts(ni);
    tl(np+1) = 10^6;
    
    % now find the t coordinate
    ip = 1;
    tx(ip) = 0.0;
    
    
    for j = 1:2 % execution loop label 5500
        xx(j,1) = x(j,1);

    end
    

    sre1 = 0.0;
    case6001 = true;

    for i = 1:nn-1 % execution loop label 6000

        sref = sre1;

        for j = 1:2 % execution loop label 6010

            x0(j) = x(j,i);
            xc(j) = c(j,i);
            x1(j) = x(j,i+1);

        end
        [srein,ax,ay] = getarcln(1.0,x0,xc,x1,srein);
        
        
        sre1 = sref + srein;
        
        
        while case6001 == true 
            
            % jump gate exit 6001
            if tl(ip+1) > sre1

                % jump gate label 6000
                break

            else

                cl = tl(ip+1) -sref;
                tint = 0.5;
                ttry = tint;

                
                for k = 1:niter % execution loop label 8000

                    tint = 0.5*tint;

                    

                    [cres, ax, ay] = getarcln(ttry,x0,xc,x1,cres);
                    % temp_x = [temp_x;t_x];
                    % temp_cres = [temp_cres; cres];

                    ctmp = cres - cl;
                    
                    

                    if ctmp <= 0

                        ttry = ttry + tint;

                    end

                    if ctmp > 0

                        ttry = ttry - tint;

                    end
                end

                ip = ip + 1;
                
                

                tx(ip) = real(i-1) + ttry;
                
                xt = getxy(ttry,x0,xc,x1);

                for j = 1:2 % execution loop label 6900

                    xx(j,ip) = xt(j);

                end


                % jump gate label 6001
            end
        end
        % jump gate exit 6000
    end

    tx(np) = real(ni-1);
    
    for j = 1:2 % execution loop label 9000
        xx(j,np) = x(j,nn);
    end
    

    return
end