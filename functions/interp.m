
function [npoin, lboud, nbno, nbou, nnn, lbou, coord, dis, alp, anx, any, t, lx,c,ti,ts,td,tx,xx, ni, np, tg, xt,arZ,i1,i2,i3,ienr,ind,tl, srein,amo,ttan, x,rcond, ax, ay, xq,ilast] = interp(npoig,neleg,npoin,coorg,ieleg,intmeg,delta,coord,coorn,nbno,nnn,lcoor,lbou,nbou,ibs,ilast,ipoii,ipoip,lboud,rcond,igraph,unit_17, t, lx,c,ti,ts,td,tx,xx, ni, np, tg,    xt,arZ,i1,i2,i3,ienr,ind,tl, srein,amo, dis,alp,anx,any,ttan, x)

    
    xc = zeros(2,1);
    
    % THIS SUBROUTINE IS CALLED ONCE FOR EACH BOUNDARY SEGMENT.
    % IT INTERPOLATES ADDITIONAL BOUNDARY POINTS ACCORDING TO THE
    % SPACING SPECIFIED BY THE BACKGROUND GRID

    % READ THE DATA 

    fid_17 = fopen(unit_17, 'r');
    fgets(fid_17);
    temp = fscanf(fid_17,' %f %f %f',[3 1]);
    temp = temp';
    nreg = single(temp(:, 1));
    nfn  = single(temp(:, 2));
    nbcs = single(temp(:, 3));
    
    fgets(fid_17);
    fgets(fid_17); % skip COORDINATES
    for i = 1:nbcs % skip COORDINATES
        fgetl(fid_17);
    end
    fgetl(fid_17); % skip BOUNDARY SEGMENTS
    
    temp = zeros(4, nbcs);
    for i = 1:nbcs
        temp(:,i) = fscanf(fid_17,'%f %f %f %f',[4 1]);
        fgets(fid_17);
        fgets(fid_17);
    end
    
    ib     = temp(1,ibs);
    nn     = temp(2,ibs);
    igeom  = temp(3,ibs);
    icond  = temp(4,ibs);

    % to single
    ib     = single(ib   );  
    nn     = single(nn   ); 
    igeom  = single(igeom); 
    icond  = single(icond); 
    coorn = single(coorn);
    % -----------

    fid_17 = fopen(unit_17, 'r');
    fgets(fid_17);
    temp= fscanf(fid_17,' %f %f %f',[3 1]);
    
    fgets(fid_17);
    fgets(fid_17); % skip COORDINATES
    for i = 1:nbcs % skip COORDINATES
        fgetl(fid_17);
    end
    fgetl(fid_17); % skip BOUNDARY SEGMENTS
    fgetl(fid_17);
    temp = zeros(2,nbcs);
    for i = 1:nbcs
        temp(:,i) = fscanf(fid_17,'%f %f',[2 1]);
        fgets(fid_17);
        fgets(fid_17);
    end
    
    ln = temp(:,ibs)';
 
    % to single
    ln = single(ln);
    % -----------
  
    % STORE THE INTERSECTION POINTS
    if (nbou == 0)

        nbou      = 2;
        npoin     = 2;
        lbou(1,1) = ln(1);
        lbou(1,2) = ln(nn);
        lbou(2,1) = 1; 
        lbou(2,2) = 2; 
        

        for j = 1:2

            rcond(j,1) = 0.0;
            rcond(j,2) = 0.0;
            coord(j,1) = coorn(j,ln(1));
            coord(j,2) = coorn(j,ln(nn));

        end
    else

        kfir = ln(1);
        ksec = ln(nn); 
        ifir = 1;
        isec = 1;

        for ibou = 1:nbou % loop 100
           
            if (lbou(1,ibou) == kfir)
                ifir = 0;
            end

            if (lbou(1,ibou) == ksec)
                isec = 0;
            end
        end

        if (ifir == 0)
            
            if (isec == 0)
                ;
            else
                nbou           = nbou + 1; 
                npoin          = npoin + 1;
                rcond(1,npoin) = 0.0;
                rcond(2,npoin) = 0.0;
                coord(1,npoin) = coorn(1,ksec);
                coord(2,npoin) = coorn(2,ksec);
                lbou(1,nbou)   = ksec;
                lbou(2,nbou)   = npoin; 
            end

        else
            nbou           = nbou + 1; 
            npoin          = npoin + 1;
            rcond(1,npoin) = 0.0;
            rcond(2,npoin) = 0.0;
            coord(1,npoin) = coorn(1,kfir);
            coord(2,npoin) = coorn(2,kfir);
            lbou(1,nbou)   = kfir;
            lbou(2,nbou)   = npoin;
        end
    end
    t = rfillm(t,2,nn,0.0);
    lx = rfilliv(lx,nn,0);

    
    % TRANSFER
    for i = 1:nn
        for j =1:2
            x(j,i) = coorn(j,ln(i));
        end
    end
   
    
    % GET THE CONTROL POINTS

    [x,c] = getcntrl(nn,x,c,lx,t,ls, xc);

    % PLOT THE CURVES
    for i = 1:nn - 1
        for j =1:2
            x0(j) = x(j  ,i);
            xc(j) = c(j  ,i);
            x1(j) = x(j,i+1);
        end

        if (igraph == 1) 
            %plotsp(x0,xc,x1)
            disp("plotsp")
        end

    end
    
    % FIND INTERSECTIONS WITH THE BACKGROUND GRID 
  

    [ti,ni, arZ, i1, i2, i3, ienr, xq,ilast] = findtin(npoig,neleg,ieleg,intmeg,coorg,ilast,nn,ti,x,c,ni,ipoii,i1,i2,i3,ienr,ind,arZ, xt);
  

    % find arc length coordinate of the intersection points

    ts = findtsn(nn,x,c,ni,ti,ts);

    % find the spacing required

  
    [td,arZ,i1,i2,i3,ienr, dis, alp, anx, any, ttan,amo,ilast]= findtdn(nn,x,c,ni,ti,td,npoig,neleg,ieleg,intmeg,coorg,delta,ilast,igeom,  amo, ttan,arZ,i1,i2,i3,ienr,dis,alp,anx,any);
        
    %find location of the boundary points

    [tx,xx,np,tl, srein, tg, x, ax, ay] = findtpo(nn,ni,np,x,c,ts,tg,td,tx,xx,tl, srein);
 
    % PLOT GENERATED POINT IF REQUIRED
    if (ipoip == 0 || igraph ~= 1)
        ;
    else
        for ip = 1:np
            ;
        end
    end

    
    nnn(ibs)     = np;
    nbno(ibs,1)  = ln(1);
    nbno(ibs,np) = ln(nn);

    % STORE ALL POINTS BUT THE END ONES
    if (np == 2) % jump gate 9000
        return
    end
    
    for ip = 2:np-1 % loop label 9001
        npoin          = npoin + 1; 
        nbno(ibs,ip)   = npoin;
        anorx          = xx(2,ip + 1) - xx(2,ip - 1);
        anory          = xx(1,ip - 1) - xx(1,ip + 1);
        amo            = sqrt(anorx^2 + anory^2);
        
        rcond(1,npoin) = anorx/amo;
        rcond(2,npoin) = anory/amo;

        for j = 1:2 % loop label 9002
            coord(j,npoin) = xx(j,ip); 
        end
        
        lboud(npoin) = icond;
        
    end
    
    return
end