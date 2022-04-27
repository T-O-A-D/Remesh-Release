function [iside,iloca,lwher,icone] = side(nelem,npoin,iloca,intmat,iside,lwher,lhowm,icone)
        
    % fill in the lhowm : nr. of elements per node


    for ip = 1:npoin % execution loop label 1490
        lhowm(ip) = 0;
    end
    
    for ie = 1:nelem % execution loop label 1500
        for in = 1:3 % execution loop label 1500
            ip = intmat(in,ie);
            lhowm(ip) = lhowm(ip) + 1;
        end
    end
    
    
    % fill in lwher : location of each node inside icone
    lwher(1) = 0;
    for ip = 2:npoin % execution loop label 1600
        lwher(ip) = lwher(ip-1) + lhowm(ip-1);
    end

    % fill in icone : elements in each node
    for ip = 1:npoin % execution loop label 1690
        lhowm(ip) = 0;
    end
    for ie = 1:nelem % execution loop label 1700
        for in = 1:3 % execution loop label 1701
            ip=intmat(in,ie);
            lhowm(ip)=lhowm(ip)+1;
            jloca=lwher(ip)+lhowm(ip);
            icone(jloca)=ie;
        end
    end
    % loop over the nodes
    iloca = 0;
    for ip = 1:npoin % execution loop label 3000
        iloc1 = iloca;
        iele = lhowm(ip);
        if iele == 0 % jump gate label 3000
            continue
        else
            % initialize iside ---> important for boundary sides
            for is = 1:iele+2 % execution loop label 3001
                iside(3,is+iloc1) = 0;
                iside(4,is+iloc1) = 0;
            end
            iwher = lwher(ip);
            % loop over the elements surrounding the point ip
            ip1 = ip;
            for iel = 1:iele % execution loop label 3090
                ie = icone(iwher + iel);
                % find out position of ip in the coneivity matrix
                for in = 1:3 % execution loop label 3091
                    in1 = in;
                    ipt = intmat(in,ie);
                    if (ipt == ip) % jump gate label 3092
                        break
                    end
                end
                % jump gate label 3092
                for j1 = 1:2 % execution label 3100
                    
                    in2 = in1 + j1;
                    if in2 > 3 
                        in2 = in2 - 3;
                    end
                    
                    ip2 = intmat(in2,ie);
                    if ip2 < ip1 % jump gate label 3100
                        ;
                    else
                        went_to_7303 = false;
                        % check the side ----> new or old
                        if iloca == iloc1 % jump gate label 7304
                            went_to_7303 = false;
                            % new side
                        else
                            for is = iloc1+1:iloca % execution loop label 5600
                                jloca = is;
                                if iside(2,is) == ip2 % jump gate label 7303
                                    went_to_7303 = true;
                                    iside(2+j1,jloca) = ie;
                                    break   
                                end
                                went_to_7303 = false;
                            end
                            
                        end 
                        % jump gate exit label 7304
                        if went_to_7303 == false
                            iloca = iloca + 1;
                            iside(1,iloca) = ip1;
                            iside(2,iloca) = ip2;
                            iside(2+j1,iloca) = ie;
                        end
                        % exit execution label 3012  
                    end

                    % jump gate exit label 3100
                    % end loop over elements surrounding point ip    
                end       
            end  

            for is = iloc1+1:iloca % execution label 8000
            
                if iside(3,is) ~= 0 % jump gate label 8000
                    continue
                else
                    iside(3,is)=iside(4,is); 
                    iside(4,is)=0;
                    iside(1,is)=iside(2,is);
                    iside(2,is)=ip1;
                end
                % jump gate exit label 8000 
            end
        % end loop over points
        end
    end

    % jump gate exit label 3000
    return
end