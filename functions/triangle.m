function [node, nelem, toler, iel,coor,fid_7 ] = triangle(nreg,nonf,npfront,nqfront,nregi,nonr,iel,coor,nelem, node,neleg,npoig,ieleg,intmeg,coorg,delta,strec,toler,ilast,igraph,unit_7,fid_7,cv,dis)
    
    iyon = single([]);
    arZ  = single([]);
    i1   = single([]);
    i2   = single([]);
    i3   = single([]);
    ienr = single([]);
    alp  = single([]);
    case603f = 1;
    text = 'aw for re-ordering, (aw=1.:always,aw=large:never)';

    if cv == '1'
        % read for aw
        fid_7 = fopen(unit_7, 'r');
        for i = 1:2
            fgets(fid_7);
        end

        temp = fscanf(fid_7,'%f',[1 1]);
        aw = temp';
    else
        aw = fscanf(fid_7,'%f',[1]);
    end    
  
   
    % ------------
    
    disp(text)
    % to single
    aw = single(aw);
    disp(aw)
    case160 = 1;
    case161 = 1;
    case360 = 1;
    a = 0;
    for in = 1:node % execution label 5

        x(in) = single(coor(1,in));
        y(in) = single(coor(2,in));

    end
    
    % loop over regions
    nelem = 0;
    dis = single(0);
    counter1 = 0;
    counter2 = 0;
    for ireg =1:nreg % execution label 100

        for ik=1:node % execution label 570

            ncheck(ik) = -100;

        end
        
        % setup ncheck values for region
        
        for ik=1:nonf(ireg) % execution label 580

            k1 = npfront(ireg,ik);
            
            ncheck(k1) = 2;

        end
        
        disw = 0.0;

        while case160 == 1 % execution loop label 160
            
            nl = single(nonf(ireg));
            counter1 = counter1 + 1;
            
            while case161 == true % execution loop label 161

                counter2 = counter2 + 1;
                kn   = single(nqfront(ireg,nl));
                kn1  = single(npfront(ireg,nl));
                xn   = single(coor(1,kn));
                yn   = single(coor(2,kn));
                xn1  = single(coor(1,kn1));
                yn1  = single(coor(2,kn1));
                alph = single(amin1([strec(2,kn1),strec(2,kn)]));
                anx  = single(0.5*(strec(3,kn1) + strec(3,kn)));
                any  = single(0.5*(strec(4,kn1) + strec(4,kn)));
                anm  = single(sqrt(anx*anx + any*any));
                anx  = anx/anm;
                any  = any/anm;
                xnf    = single(fxba(xn,yn,anx,any,alph));
                ynf    = single(fyba(xn,yn,anx,any,alph));
                xn1f   = single(fxba(xn1,yn1,anx,any,alph));
                yn1f   = single(fyba(xn1,yn1,anx,any,alph));
                x12f   = single(xnf - xn1f);
                y12f   = single(ynf - yn1f);
                alen1f = single(sqrt(x12f*x12f + y12f*y12f));
                
                if (alen1f < aw*disw)
                    ; %% goto 162
                    
                else

                    disw1 = 0.0;
                    [npfront, nqfront, disw] = order(nl,ireg,npfront,nqfront,coor,npoig,neleg,coorg,ieleg,delta,strec,disw,disw1);
                    %npfront
                    %nqfront
     
                    continue %GOTO 161

                end
                
                
                %come to 162
                tole1 = single(0.00001*alen1f);
                
                % find out average spacing
                average = amin1([strec(1,kn1),strec(1,kn)]);
                
                x12     = single(xn - xn1);
                y12     = single(yn - yn1);
                alen1   = single(sqrt(x12*x12 + y12*y12));
    
               
                % create a new node
                % safety factor
                csafe = single(1);
                dside = single(csafe*average);

                if dside > single(2.0*alen1f)
                    dside  = single(2.0*alen1f);
                end

                if dside < single(0.55*alen1f)
                    dside = single(0.55*alen1f);
                end

                twod2  = single(2.0*dside*dside);
                xbarf  = single(0.5*(xnf+xn1f));
                ybarf  = single(0.5*(ynf+yn1f));
                xdiff  = single(xnf-xn1f);
                ydiff  = single(ynf-yn1f);
                dkarg  = single(xdiff*xdiff+ydiff*ydiff);
                d12    = single(sqrt(dkarg));
                hkarg  = single(0.5*twod2-0.25*d12*d12);
                hk     = single(sqrt(hkarg));
                xcbf   = single(-hk*ydiff/d12);
                ycbf   = single(hk*xdiff/d12);
                xtempf = single(xbarf+xcbf);
                ytempf = single(ybarf+ycbf);
                
                %xtemp=xtr(xtempf,ytempf,anx,any,alph);
                %ytemp=ytr(xtempf,ytempf,anx,any,alph);

                xtemp  = fxtr(xtempf,ytempf,anx,any,alph);
                ytemp  = fytr(xtempf,ytempf,anx,any,alph);
                
                %loop over possible nodes - find closest neighbours
                a      = single(yn1-yn);
                b      = single(xn-xn1);
                c      = single((xn1-xn)*yn1+(yn-yn1)*xn1);
                radius = single(dside);
                h1     = single(0.8*radius);
                inum   = single(0);

                
                % loop 110
                for kp = 1:nonr(ireg)

                    ken = single(nregi(ireg,kp));
                    
                    if single(ken == kn) || single(ken == kn1)
                        continue
                        %case110 = 1;
                    else

                        %case110 = 0;

                        xken=single(coor(1,ken));
                        yken=single(coor(2,ken));
                        %xkenf=xba(xken,yken,anx,any,alph);
                        %ykenf=yba(xken,yken,anx,any,alph);
                        xkenf=single(fxba(xken,yken,anx,any,alph));
                        ykenf=single(fyba(xken,yken,anx,any,alph));
                        xdiff1=single(xkenf-xtempf);
                        ydiff1=single(ykenf-ytempf);
                        distf=single(sqrt(xdiff1*xdiff1+ydiff1*ydiff1));

                        if single(distf > h1)
                            continue
                            %case110 = 1;
                        else
                            %case110 = 0;
                            if single(((a*xken) + (b*yken) + c)) <= 0.0
                                %case110 = 1;
                                continue
                            else

                                %case110 = 0;
                                inum=single(inum+1);
                                howf1(inum)=single(distf);
                                near1(inum)=single(ken);

                            end
                        end
                    end
                end
                %end

                % decide which of the nodes is chosen
                if inum == 0

                    inum = 1;
                    near(1) = 0;
                    howf(1) = 0.0;

                else
                    % order them
                    for i = 1:inum % execution 599

                        comp = single(1*10^6);
                        for j = 1:inum % execution 598

                            if near1(j) == 0
                                continue
                                % GOTO 598
                            else
                                if howf1(j) > comp
                                    continue
                                    % GOTO 598
                                else
                                    is = j;
                                    comp=howf1(j);
                                end
                            end
                        end

                        near(i) = single(near1(is));
                        howf(i) = single(howf1(is));
                        near1(is) = 0;
                        
                    end

                    %add the new point to the list
                    inum = single(inum + 1);
                    near(inum) = 0;
                    howf(inum) = 0.0;

                end
                % select---> start by the closest
                %case601 = 1;
                %while case601 == 1
                for i = 1:inum % execution label 601

                    kp=near(i);

                    if kp == 0

                        xp=single(xtemp);
                        yp=single(ytemp);

                    else

                        xp=single(coor(1,kp));
                        yp=single(coor(2,kp));

                    end

                    % see if this connection is possible
                    
                    iyon = possible(kn1,kn,kp,xn1,yn1,xn,yn,xp,yp,ireg,nonf,nonr,npfront,nqfront,nregi,coor,iyon);
                    
                    if iyon == 0
                        
                        case603 = 0;
                        continue
                    else
                        case603 = 1;
                        break
                    end

                end

                % we are in trouble !!!!!
                % find the 30 existing nodes that give maximum angle
                if case603 == 0
                    
                    inum = 0;
                    ang1 = 0.0;
                    for kp = 1:nonr(ireg) % execution label 210
                        ken = single(nregi(ireg,kp));
                        if ken == kn || ken == kn1 % jump gate label 210
                            continue
                        else
                            xken = single(coor(1,ken));
                            yken = single(coor(2,ken));
                            if single((a*xken) + (b*yken) + c) <= 0.0 % jump gate label 210
                                continue
                            else
                                % see if this connection is possible

                                iyon = possible(kn1,kn,ken,xn1,yn1,xn,yn,xken,yken,ireg,nonf,nonr,npfront,nqfront,nregi,coor,iyon);
                                if (iyon == 0) % jump gate label 210
                                    
                                    continue
                                else

                                    xkenf  = single(fxba(xken,yken,anx,any,alph));
                                    ykenf  = single(fyba(xken,yken,anx,any,alph));
                                    xdiff1 = single(xkenf-xn1f);
                                    ydiff1 = single(ykenf-yn1f);
                                    xdiff2 = single(xkenf-xnf);
                                    ydiff2 = single(ykenf-ynf);
                                    distf1 = single(sqrt(xdiff1*xdiff1+ydiff1*ydiff1));
                                    distf2 = single(sqrt(xdiff2*xdiff2+ydiff2*ydiff2));
                                    cosa   = single((xdiff1*xdiff2+ydiff1*ydiff2)/(distf1*distf2));

                                    if (cosa > 1.0)
                                        cosa = single(1.0);
                                    end
                                    
                                    if cosa < -1.0
                                        cosa = single(-1.0);
                                    end

                                    angl = single(acos(cosa));
                                    
                                    if angl < ang1 % jump gate label 210
                                        
                                        continue

                                    else
                                        
                                        howf(inum + 1) = single(angl);
                                        near(inum + 1) = single(ken);

                                        if inum == 0 % jump gate label 311
                                            ;
                                        else
                                            
                                            for i = 1:inum % execution label 310

                                                k = i;
                                                angi = single(howf(i));

                                                if angi > angl

                                                    continue % jump gate label 310
                                                else
                                                    
                                                    for j = 1:(inum-k) % execution label 410
                                                        l = inum - j;
                                                        near(l+1) = single(near(l));
                                                        howf(l+1) = single(howf(l));
                                                    end

                                                    near(k) = single(ken);
                                                    howf(k) = single(angl);
                                                    
                                                    break
                                                    % jump gate label 310
                                                end
                                            end
                                            % jump gate exit label 310
                                            
                                            
                                        end
                                        if (inum < 30)

                                            inum = inum + 1;
                                            
                                        end

                                        if (inum == 30)

                                            ang1 = single(howf(30));

                                        end
                                    end
                                end
                            end
                        end
                    end

                    % selection ---> start by closest
                    
                    if inum == 0
                        case603f = true;
                        % jump gate label 703

                    else

                        wfar = single(1*10^6);
                        for i = 1:inum % execution label 701

                            kp = single(near(i));
                            angi = single(howf(i));
                            xp = single(coor(1,kp));
                            yp = single(coor(2,kp));

                            % see if this connection is possible
                            % call possible(kn1,kn,kp,xn1,yn1,xn,yn,xp,yp,ireg,nonf,
                            % nonr,npfront,nqfront,nregi,coor,iyon)
                            % if(iyon.eq.0) goto 701
                            % if(angi.gt.0.7853981) goto 603
                            % check now the size


                            xpf = fxba(xp,yp,anx,any,alph);
                            ypf = fyba(xp,yp,anx,any,alph);
                            d1  = single(sqrt(((xpf-xn1f)^2)+((ypf-yn1f)^2)));
                            d2  = single(sqrt(((xpf-xnf)^2)+((ypf-ynf)^2)));
                            di  = single(amax1([d1,d2]));

                            if (di < wfar)
                                kp1 = single(kp);
                            end
                            
                            if (di < wfar)
                                wfar = single(di);
                            end
                        end
                        
                        kp = single(kp1);
                        xp = single(coor(1,kp));
                        yp = single(coor(2,kp));
                        
                        case603f = false; % jump gate label 603f              
                    end

                    % jump gate exit label 703
                    if case603f == true
                        msg = 'cannot find the connectivity';
                        fprintf(msg);
                        stop
                        % try another side
                        disw1 = 1.01*alen1f;

                        [npfront, nqfront, disw] = order(nl,ireg,npfront,nqfront,coor,npoig,neleg,coorg,ieleg,delta,strec,disw,disw1);
                        case161 = 1;
                        
                        continue; % GOTO 603
                            
                    end
                end
                % nothing wrong with it 
                % jump gate label 603
                indic = 1;
                knear = single(kp);

                if knear ~= 0 
                    indic = 0;
                end

                if indic == 0 % jump gate label 620

                else

                    node = node + 1;
                    coor(1,node) = single(xp);
                    coor(2,node) = single(yp);

                    % find out interpolated variables
                    inorm = 1;

                    [ienr, arZ, i1, i2, i3] = findel(npoig,neleg,coorg,ieleg,intmeg,xp,yp,ilast,arZ,i1,i2,i3,ienr);

                    [dis, alp, anx, any] = getvalue(npoig,i1,i2,i3,arZ,delta,dis,alp,anx,any,inorm);

                    strec(1,node) = single(dis);
                    strec(2,node) = single(alp);
                    strec(3,node) = single(anx);
                    strec(4,node) = single(any);

                    nonr(ireg) = nonr(ireg) + 1;
                    nuno = nonr(ireg);
                    nregi(ireg,nuno) = node;
                    knear = node;
                    ncheck(knear) = -100;
                end

                % jump gate label 620

                % form element
                nelem = nelem+1; 
                iel(1,nelem) = single(kn1);
                iel(2,nelem) = single(kn);
                iel(3,nelem) = single(knear);

                % UPDATE FRONT AND ACTIVE NODES
                nl2=single(nl+1);
                nqfront(ireg,nl2) = single(kn);
                nqfront(ireg,nl)  = single(knear);
                npfront(ireg,nl2) = single(knear); 

                if ncheck(knear) < 0
                    ncheck(knear) = 0;
                end

                ncheck(knear) = ncheck(knear) + 2;
                nonf(ireg) = nl2;

                % delete sides from active list
                ncht   = 0;
                npass  = 1;
                ntop   = single(kn1);
                nbot   = single(knear);
                marker = 0;
                knon   = single(nonf(ireg) - 2);

                % come to 360
                while case360 == true

                    for kp = 1:knon % execution label 300

                        ktop = single(npfront(ireg,kp));
                        kbot = single(nqfront(ireg,kp));

                        if marker > 0 

                            kp1 = single(kp - 1);
                            npfront(ireg,kp1) = single(ktop);
                            nqfront(ireg,kp1) = single(kbot);

                        else

                            if (ktop == nbot && kbot == ntop) 

                                marker = 1;
                                nonf(ireg) = single(nonf(ireg) - 2);
                                ncheck(ktop) = single(ncheck(ktop) - 2);
                                ncheck(kbot) = single(ncheck(kbot) - 2);
                                ncht = 1;

                            else
                                
                                continue % jump gate label 300

                            end
                        end
                    end

                    npass = npass + 1;

                    if npass > 2 % jump gate label 340
                        % jump gate label 340
                        break;       
                    else

                        if marker == 0 % jump gate label 350

                            knon = knon + 1; % jump gate label 350

                        else

                            npfront(ireg,nonf(ireg)) = single(knear);
                            nqfront(ireg,nonf(ireg)) = single(kn);
                            marker = 0;
                            knon = knon - 1;

                        end
                        % label 400
                        ntop = knear;
                        nbot = kn;
                        continue; % GOTO 360
                    end
                end

                % remove nodes from active list
                if ncht == 0 % jump gate label 370

                else

                    ired = 0;
                    knon = nonr(ireg);

                    for kp = 1:knon % execution loop label 380

                        kkn = single(nregi(ireg,kp));
                        kch = single(ncheck(kkn));

                        if (kch ~= 0) % jump gate label 390

                        else

                            ired = ired + 1;
                            continue %GOTO 380

                        end

                        kp1 = kp - ired; % jump gate label 390
                        nregi(ireg,kp1) = nregi(ireg,kp);

                    end

                    nonr(ireg) = knon - ired;
                end

                % jump gate label 370
                if (nonf(ireg) > 0) % jump gate label 160   
                    case160 = true;
                    break; %GOTO 160

                else
                    case160 = 0;
                    break; %GOTO 100
                end
                
                % end of loop over regions
            
            end

            if case160 == 0
                break; %GOTO 100
            end

        end

    end

    % jump gate exit label 100

    return
end
