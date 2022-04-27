function [ti, ni, arZ, i1, i2, i3, ienr, xq,ilast] = findtin(npoig,neleg,ieleg,intmeg,coorg,ilast, nn,ti,x,c,ni,ipoii,i1,i2,i3,ienr,ind,arZ, xt)

    %THIS SUBROUTINE FINDS OUT THE INTERSECTION POINTS OF A BEZIER
    %QUADRATIC CURVE WITH THE BACKGROUND GRID MESH
    
    pass1000 = false;
    % FIRST POINT
    t = single(0.0);
    ni = single(1);
    ti(ni) = single(0.0);
    % goto 2020 help bool
    % LOOP OVER THE NUMBER OF LOCAL SEGMENTS
    % loop 1000

    for i = 1:nn-1 
        tref = real(i-1);
        
        % loop 1001
        for j = 1:2
            x0(j)=x(j,i);
            xc(j)=c(j,i);
            x1(j)=x(j,i+1);
        end

        
        % loop 2020
        while true

            t=amax1([(ti(ni)-tref+0.0001),0.0]);
            xt = getxy(t,x0,xc,x1);
            
            [ienr, arZ, i1, i2, i3,ilast] = findel(npoig,neleg,coorg,ieleg,intmeg,xt(1),xt(2),ilast,arZ,i1,i2,i3,ienr);

            % loop 1004
            for j = 1:3
                j1=j+1;
                j2=j+2;

                if (j1 > 3)
                    j1=j1-3;
                end

                if (j2 > 3)
                    j2=j2-3;
                end
                ip=ieleg(j1,ienr); 
                iq=ieleg(j2,ienr);
                
                % loop 1005
                for k = 1:2
                    xp(k)=coorg(k,ip);
                    xq(k)=coorg(k,iq);
                end

                tmin=amax1([(ti(ni)-tref+0.0001),0.0]);
                
                [ind, t] = cross(x0,xc,x1,xp,xq,tmin,ind,t);

                if (ind == 0)
                    tpos(j)=1.e+6;
                else
                    tpos(j)=tref+t;
                end
            end

            % ANY POINTS ?
            tpmin = amin1([tpos(1),tpos(2),tpos(3)]);

            if (tpmin > 1.e+5)
                
                % goto 990
                break
            else
                % YES !!!
                ni=ni+1;
                ti(ni)=tpmin;
                xt = getxy(tpmin-tref,x0,xc,x1);

                % loop 1100
                for j = 1:3
                    
                    if (tpos(j) < 1.01*tpmin)
                        k = j;
                    end
                end

                ilast = intmeg(k,ienr);

                % goto 2020
            end
        end

        % 990 CONTINUE
        % ADD T=0.5 IF NECESSARY

        tcomp=tref+0.5;
        
        % loop 970
        for ii = 1:ni
            
            if (abs(ti(ii)-tcomp) < 0.3)
                pass1000 = true;

                % goto 1000
                break
            else
                pass1000 = false;
            end
        end

        if pass1000 == true

            continue

        end

        % 970 CONTINUE
        ni=ni+1;
        ti(ni)=tcomp;

        % loop 980
        for ii = 1:ni
            
            if (ti(ni-ii) < tcomp)
                % goto 981
                break
            end

            ti(ni-ii+1)=ti(ni-ii);
            ti(ni-ii)=tcomp;
        end

        % 981 CONTINUE
        xt = getxy(0.5,x0,xc,x1);
    end
    % 1000 CONTINUE

    % ADD THE LAST POINT

    if (ti(ni) > tref+0.999)
        return
    end

    ni=ni+1;
    ti(ni)=tref+1;

    return
end