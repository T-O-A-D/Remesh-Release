function [td,arZ,i1,i2,i3,ienr, dis, alp, anx, any, ttan,amo,ilast] = findtdn(nn,x,c,ni,ti,td,npoig,neleg,ieleg,intmeg,coorg,delta,ilast,igeom,  amo, ttan,arZ,i1,i2,i3,ienr,dis,alp,anx,any) 

    atole = single(0.1);
    % loop over number of points
    for ii = 1:ni % execution loop label 1000

        tot = ti(ii) - 0.00001;
        in  = fix(tot);   
        t   = tot - real(in) + 0.00001;
        in  = in + 1;
        
        % transfer
        for j = 1:2 % execution loop label 1001

            x0(j) = x(j,in); 
            xc(j) = c(j,in); 
            x1(j) = x(j,in+1);

        end
        
        % get the coordinates
        xt = getxy(t,x0,xc,x1);
        % get the tangent direction
        [ttan,amo] = gettan(t,x0,xc,x1,ttan,amo);

        % the curvature
        curv = 0.0;
        if igeom == 1 
            curv = getcur(t,x0,xc,x1,curv); 
        end

        % now interpolate the spacing
        inorm = 1;
        [ienr, arZ, i1,i2,i3,ilast] = findel(npoig,neleg,coorg,ieleg,intmeg,xt(1),xt(2),ilast,arZ,i1,i2,i3,ienr);
    
        [dis, alp, anx, any] = getvalue(npoig,i1,i2,i3,arZ,delta,dis,alp,anx,any,inorm);


        % fill in id

        xfi = fxba(ttan(1),ttan(2),anx,any,alp);
        yfi = fyba(ttan(1),ttan(2),anx,any,alp);
        afi = sqrt(xfi*xfi + yfi*yfi);
        xfi = xfi/afi;
        yfi = yfi/afi;
        xre = fxtr(xfi,yfi,anx,any,alp);
        yre = fytr(xfi,yfi,anx,any,alp);
        red = sqrt(xre*xre+yre*yre);
        dis = dis*red;
        
        if igeom ~= 1 % jump gate label 2000

            % goto 2000

        else

            if curv*dis < atole % jump gate label 2000

                % goto 2000

            else

                dis = atole/curv;

            end
        end
        % jump gate exit label 2000
        td(ii) = 1./dis;
    end

    return
end

