function mmat = gtmmat(nelem,npoin,nnode,ngeom,intmat,geome,mmat)
    c6 = 1/6;
    % zero mmat
    mmat = rfillv(mmat,npoin,0.0);
    % loop over the elements
    for ielem = 1:nelem % execution loop label 1000
        % jacobian of the element
        rj = geome(7,ielem);
        rj6 = rj*c6;
        % add to mmat
        for inode = 1:nnode % execution loop label 1020
            in = intmat(inode,ielem);
            mmat(in) = mmat(in) + rj6;
        end
        % end of loop over the elements
    end
    % control output of mmat
    % inversion of mmat
    for i = 1:npoin % execution loop label 2000
        mmat(i) = 1/mmat(i);
    end
    
    return
end