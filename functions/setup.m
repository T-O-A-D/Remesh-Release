function [nonf,nonr,nregi,npfront,nqfront] = setup(nreg,mbcs,ibcs,nbno,npoin,coor,node,npfront,nqfront,nonf,nregi,nonr)

    % SET UP GENERATION FRONTS AND NODES - REGION BY REGION


    % LOOP OVER REGIONS

    % loop 100
    for ireg = 1:nreg

        % COLLECT REGION BOUNDARY NODES
        % LOOP OVER BOUNDARY SEGMENTS 
        lnode=0;

        % loop 520
        for iseg = 1:mbcs(ireg)

            noseg = abs(ibcs(ireg,iseg));
            np1   = nbno(noseg,1);
            nn    = npoin(noseg);
            np2   = nbno(noseg,nn);
            nl1   = nn - 1;

            if (ibcs(ireg,iseg) < 0)

                % goto 540

                % loop 570
                for kn = 1:nl1

                    in    = nl1 + 1 - kn;
                    lnode = lnode + 1; 
                    npfront(ireg,lnode) = nbno(noseg,in + 1);
                    nqfront(ireg,lnode) = nbno(noseg,in);
                    
                end

            else

                % loop 550
                for kn = 1:nl1

                    lnode = lnode + 1; 
                    npfront(ireg,lnode) = nbno(noseg,kn);
                    nqfront(ireg,lnode) = nbno(noseg,kn + 1);

                end

                % goto 560
            
            end

        % END OF LOOP OVER BOUNDARY SEGMENTS
        end

        % SET UP INITIAL SET OF NODES FOR REGION

        % loop 310
        for knode = 1:lnode
            nregi(ireg,knode) = npfront(ireg,knode);
        end

        nonr(ireg) = lnode;
        nonf(ireg) = lnode;


    %END LOOP OVER REGIONS
    end

    return
end