function [iside, lcore, intmat,iloca] = swapdi(npoin,nelem,nside,iside,intmat,lcore,lcoid,lcoor,coord,lwher,lhowm,icone)
    iswaplp = 0;
    i3 = int16(32767);
    while true
        ichan = single(0);
        % loop over the sides 
        for is = 1:nside % execution loop label 2000
            i1  = int16(iside(1,is));
            i2  = int16(iside(2,is));
            ie1 = int16(iside(3,is));
            ie2 = int16(iside(4,is));
            % check for boundary sides
            if ie2 == 0 % jump gate label 2000

                continue

            else
                % determine i3 & i4
                for in = 1:3 % execution loop label 5000
                    in1 = in + 1;
                    if in1 >3   
                        in1 = in1 - 3;
                    end
                        in2 = in + 2;
                    if in2 > 3
                        in2 = in2 - 3;
                    end
                    ip11 = intmat(in,ie1);
                    ip12 = intmat(in1,ie1);
                    ip13 = intmat(in2,ie1);
                    if (ip11 == i1 && ip12 == i2)
                        i3 = ip13;
                    end
                    ip21 = intmat(in,ie2);
                    ip22 = intmat(in1,ie2);
                    ip23 = intmat(in2,ie2);
                    if (ip21 == i2 && ip22 == i1)
                        i4 = ip23;
                    end
                end

                % find deficit or superhavit of connectivities
                ih1  = abs(lcore(i1) - lcoid(i1));

                ih2  = abs(lcore(i2) - lcoid(i2));

                ih3  = abs(lcore(i3) - lcoid(i3));

                ih4  = abs(lcore(i4) - lcoid(i4));

                ihf1 = abs(lcore(i1) - 1 - lcoid(i1));

                ihf2 = abs(lcore(i2) - 1 - lcoid(i2));

                ihf3 = abs(lcore(i3) + 1 - lcoid(i3));

                ihf4 = abs(lcore(i4) + 1 - lcoid(i4));

                % check if it is worth to swap
                iswap = 0;
                iactu = ih1  + ih2  + ih3  + ih4;
                ifutu = ihf1 + ihf2 + ihf3 + ihf4;

                if iactu > ifutu 

                    iswap = 1;

                end

                iam = max([ih1 ih2 ih3 ih4]);
                ifm = max([ihf1 ihf2 ihf3 ihf4]);

                % disp(ih1)
                if ((iactu == ifutu) && (iam > (ifm +1))) 

                    iswap = 1;

                end

                if iswap == 0
                    %disp("iswap == 0")
                    continue
                    % jump gate label 2000
                else

                    %check area
                    x1=coord(1,i1);
                    x2=coord(1,i4);
                    x3=coord(1,i3);
                    y1=coord(2,i1);
                    y2=coord(2,i4);
                    y3=coord(2,i3);
                    %if(deter(x1,y1,x2,y2,x3,y3).le.0.0) goto 2000
                    dtmp = deter(x1,y1,x2,y2,x3,y3);
                    if dtmp <= 0.0 

                        disp("upper")
                        continue
                    % jump gate label 2000
                    else
                        x1=coord(1,i3);
                        x2=coord(1,i4);
                        x3=coord(1,i2);
                        y1=coord(2,i3);
                        y2=coord(2,i4);
                        y3=coord(2,i2);
                        dtmp = deter(x1,y1,x2,y2,x3,y3);
                        if dtmp <= 0 
                            disp("lower")
                            continue
                            % jump gate label 2000
                        else
                            % swap
                            ichan=ichan+1; 
                            intmat(1,ie1)=i1;
                            intmat(2,ie1)=i4;
                            intmat(3,ie1)=i3;
                            intmat(1,ie2)=i3;
                            intmat(2,ie2)=i4;
                            intmat(3,ie2)=i2;
                            lcore(i1)=lcore(i1)-1;
                            lcore(i2)=lcore(i2)-1;
                            lcore(i3)=lcore(i3)+1;
                            lcore(i4)=lcore(i4)+1;
                            %detect the sides i1-i4 and i3-i2
                            ist1 = 0;
                            ist2 = 0;
                            for ist = 1:nside % execution loop label 4000
                                i1t = iside(1,ist);
                                i2t = iside(2,ist);

                                if (i1t == i1 && i2t == i4) 

                                    ist1 = ist;

                                end
                                if (i1t == i4 && i2t == i1)

                                    ist1 = ist;

                                end
                                if (i1t == i3 && i2t == i2)

                                    ist2 = ist;

                                end

                                if (i1t == i2 && i2t == i3)

                                    ist2 = ist;

                                end

                            end % execution loop label 4000 end

                            if (ist1 ~= 0 && ist2 ~= 0) % jump gate label 4001
                                % jump gate label 4001
                            else
                                msg = 'error in swapdi\n';
                                %disp(ist)
                                fprintf(msg)
                                return
                                % stop function
                            end

                            % jump gate exit label 4001
                            if(iside(3,ist1) == ie2) 

                                iside(3,ist1)=ie1;

                            end
                            if(iside(4,ist1) == ie2) 

                                iside(4,ist1)=ie1;

                            end
                            if(iside(3,ist2) == ie1) 

                                iside(3,ist2)=ie2;

                            end
                            if(iside(4,ist2) == ie1) 

                                iside(4,ist2)=ie2;

                            end
                            % update iside
                            iside(1,is) = i4;
                            iside(2,is) = i3;
                        end 
                    end
                end
            % jump gate exit label 2000
            end
        end
        fprintf('%3f sides have been swapped\n',ichan)
        iswaplp = iswaplp + 1;

        if ichan >= 1 % jump gate label 1000
            %fill in iside   
            continue
        else
            break
        end
    end
        if iswaplp > 1 
            msg = 'filling iside\n';
            fprintf(msg);
            [iside,iloca,lwher,icone] = side(nelem,npoin,nside,intmat,iside,lwher,lhowm,icone);
        else
            iloca = nside;
        end
        
    return
end