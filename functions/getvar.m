function [deri2 nx ny fid_7] = getvar(nnode,namat,ngeom,nderv,nelem,npoin,coord,intmat,derip,deri2,deria,geome,mmat,unkno,ivari,fid_7)
    %initialize
    derip = rfillm(derip,nderv,npoin,0.0);
    deri2 = rfillm(deri2,4,npoin,0.0);
    %obtain mmat (already inverted)
    mmat = gtmmat(nelem,npoin,nnode,ngeom,intmat,geome,mmat);
    unkno = single(unkno);

    %read gamma
    if ivari < 3
        % jump gate label 22
        gamma = 0;
    else
        gamma = input(' enter the value of gamma (usually 1.4)');
    end

    % jump gate exit label 22
    % loop over the elements
    for ielem = 1:nelem % execution loop label 7000
        roelx = 0.0;
        roely = 0.0;
        % obtain the values needed
        for inode = 1:nnode % execution loop label 7100
            ipoin =intmat(inode,ielem);
            ro = unkno(1,ipoin);
            u1 = unkno(2,ipoin);
            u2 = unkno(3,ipoin);
            en = unkno(4,ipoin);
            vsq = u1*u1+u2*u2;
            pre = (gamma-1.)*ro*(en-0.5*vsq);
            if ivari == 1
                ropoi = ro;
            end
            if ivari == 2
                ropoi = sqrt(vsq);
            end
            if ivari == 3
                ropoi = pre/(ro^gamma);
            end
            if ivari == 4
                ropoi = sqrt(vsq*ro/(gamma*pre));
            end
            rnx = geome(inode,ielem);
            rny = geome(inode+nnode,ielem);
            roelx = roelx+rnx*ropoi;
            roely = roely+rny*ropoi;
            node(inode) = ipoin;
            nx(inode) = rnx;
            ny(inode) = rny;
        end
        % element area
        rj = geome(ngeom,ielem)*0.5/3.0;
        % assembly 
        for inode = 1:nnode % execution loop laabel 1101
            ipoin = node(inode);
            derip(1,ipoin) = derip(1,ipoin)+rj*roelx;
            derip(2,ipoin) = derip(2,ipoin)+rj*roely;
        end
        % end of loop over the elements
    end
    % multiply by mmat(inverted and divide by the nr.of el.per node
    for ipoin = 1:npoin
        for j = 1:2
            derip(j,ipoin) = derip(j,ipoin)*mmat(ipoin);
        end
    end
    % obtain the second variations
    % loop over the elements
    for ielem = 1:nelem % execution loop label 8000
        roexx = 0.0;
        roexy = 0.0;
        roeyx = 0.0;
        roeyy = 0.0;
        % obtain the values needed
        for inode = 1:nnode % execution loop label 8100
            ipoin = intmat(inode,ielem);
            ropox = derip(1,ipoin);
            ropoy = derip(2,ipoin);
            rnx   = geome(inode      ,ielem);
            rny   = geome(inode+nnode,ielem);
            roexx = roexx+rnx*ropox;
            roexy = roexy+rny*ropox;
            roeyx = roeyx+rnx*ropoy;
            roeyy = roeyy+rny*ropoy;
            node(inode) = ipoin;
            nx(inode) = rnx;
            ny(inode) = rny;
        end
        % element area
        rj = geome(ngeom,ielem)*0.5/3.0;
        % assembly
        for inode = 1:nnode % execution loop label 8101
            ipoin=node(inode);
            deri2(1,ipoin)=deri2(1,ipoin)+rj*roexx;
            deri2(2,ipoin)=deri2(2,ipoin)+rj*roexy;
            deri2(3,ipoin)=deri2(3,ipoin)+rj*roeyx;
            deri2(4,ipoin)=deri2(4,ipoin)+rj*roeyy;
        end
        % end of loop over the elements
    end
    % multiply by mmat(inverted) and divide by the nr.of el.per node
    for ipoin = 1:npoin % execution loop label 8001
        for j = 1:4 % execution loop label 8002
            deri2(j,ipoin)=deri2(j,ipoin)*mmat(ipoin);
        end
    end
    % do some smoothing if required
    %nsmoo = input('How many smoothing loops?');
    disp('How many smoothing loops?')
    nsmoo = fscanf(fid_7,'%f',[1 1]);
    if nsmoo == 0 
        % jump gate label 7891
    else
        for is = 1:nsmoo % execution loop label 7890
            deria = rfillm(deria,4,npoin,0.0);
            for ielem = 1:nelem % execution loop label 7900
                rj=geome(ngeom,ielem)*0.5/3;
                r1=0.0;
                r2=0.0;
                r3=0.0;
                r4=0.0;
                for inode = 1:nnode % execution loop label 7901
                    ipoin=intmat(inode,ielem);
                    node(inode)=ipoin;
                    r1=r1+deri2(1,ipoin)*rj;
                    r2=r2+deri2(2,ipoin)*rj;
                    r3=r3+deri2(3,ipoin)*rj;
                    r4=r4+deri2(4,ipoin)*rj;
                end
                % distribute
                for inode = 1:nnode % execution loop label 7902
                    ipoin = node(inode);
                    deria(1,ipoin)=deria(1,ipoin)+r1; 
                    deria(2,ipoin)=deria(2,ipoin)+r2; 
                    deria(3,ipoin)=deria(3,ipoin)+r3; 
                    deria(4,ipoin)=deria(4,ipoin)+r4; 
                end
            end
            % multiply by mmat(invertd) and divide by the nr. of el.per node
            for ipoin = 1:npoin % execution loop label 9001
                for j = 1:4 % execution loop label 9002
                    deri2(j,ipoin) = deria(j,ipoin)*mmat(ipoin);
                end
            end
        end
    end
    % symmetrize by equate cross derivatives
    for ipoin = 1:npoin % execution loop label 8060
        tauxy=(deri2(2,ipoin)+deri2(3,ipoin))/2;
        deri2(2,ipoin)=tauxy;
        deri2(3,ipoin)=tauxy;
    end
    
    return
end
            
        