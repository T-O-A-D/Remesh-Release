function geome = getgeo(nelem,npoin,nnode,ngeom,intmat,coord,geome,nstop)

    % this function evaluates n,x & n,y for each element shape function (linear
    % triangles) and the jacobian (=2*area)
    % loop over the elements
    nxi = [-1.0 1.0 0.0];
    net = [-1.0 0 1.0];
    for ielem = 1:nelem % execution loop label 1000
        % pick up the values needed
        for inode = 1:nnode % execution loop label 1001
            in = intmat(inode,ielem);
            x(inode) = coord(1,in);
            y(inode) = coord(2,in);
        end
        % evaluate the geometrical quantities needed
        x21 = x(2)-x(1);
        x31 = x(3)-x(1);
        y21 = y(2)-y(1);
        y31 = y(3)-y(1);
        rj  = x21*y31-x31*y21;
        if rj > 0.0 % jump gate label 2030
        else
            % print 12, i elem
            fprinf(' Some element have got negative area')
            nstop = 1;
            continue
        end
        % jump gate exit label 2030
        rj1 = 1/rj;
        xix = y31*rj1; 
        xiy =-x31*rj1; 
        etx =-y21*rj1; 
        ety = x21*rj1; 
        %----form n,x & n,y
        for in = 1:3 % execution loop label 1002
            rnxi = nxi(in); 
            rnet = net(in); 
            geome(in  ,ielem) = xix*rnxi+etx*rnet;
            geome(in+3,ielem) = xiy*rnxi+ety*rnet;
        end
    geome(7,ielem) = rj;
    % end of loop over the elements
    end
    
    return
end