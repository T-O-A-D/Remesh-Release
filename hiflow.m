%  Program hiflow
%
%  A cell-centered computer program for solving the Euler equations
%  for two-dimensional inviscid high-speed compressible flow.
%                          Professor Dr. Pramote Dechaumphai
%                          Witthaya Sadubsarn
%                          THANAPOL SUPITAYAKUL
%                          Faculty of Engineering
%                          Chulalongkorn University

%
clear 
close all
warning off;
addpath('functions')

while (true)
    inp = input('.(d)at or .(f)ix data input?: ', 's');
    
    % handle wrongly typed input
    if inp == 'd' || inp == 'D' || inp == 'f' || inp == 'F'     
        break
    else
        continue;
    end
end
cv = 0;
if inp == 'f'
        filenameGEO = input('Enter the input geometry file name: ', 's');
        cv = 1;
        % Generate mesh from background grid
        Generate(filenameGEO);
        while (true)
            cont = input('Continue:? (Y/N) ', 's');
            
            % handle wrongly typed input
            if cont(1) == 'y' || cont(1) == 'Y'        
                break
            elseif cont(1) == 'n' || cont(1) == 'N'        
                return
            else
                continue;
            end
        end
        
        filenameDAT = filenameGEO;
        cv_s = int2str(cv);
else
        filenameDAT = input('Enter the input data file name: ', 's');
        filenameGEO = filenameDAT;
        cv = input('Enter file version:');
        cv_s = int2str(cv);
end
clf('reset')
%filename = input('Enter the input file name: ', 's');
%filename = append(filename,'.dat');
%filename = "2 Shock over wedge.dat";
%filename = "7 Flow over diamond shape surface.dat";
csafe = 0.5;
while 1
    try
        
        fclose('all');
        filename = append(filenameDAT,'.dat',  cv_s);
        %filename = append(filenameGEO,'.dat');
        fid = fopen(filename, 'r');

        %
        %  Read the title of computation:
        %

        clf('reset')
        ntitle = fscanf(fid,'%i',1);
        for i = 1:ntitle+1
            text1 = fgetl(fid);
        end
        %

        %  Read input data:
        %
        text2 = fgetl(fid);
        nelem = fscanf(fid,'%i',1);
        npoin = fscanf(fid,'%i',1);
        nboun = fscanf(fid,'%i',1);
        niter = fscanf(fid,'%i',1);
        nshow = fscanf(fid,'%i',1);
        %

        for i=1:2
            text3 = fgetl(fid);
        end
        gamma = fscanf(fid,'%f',1);
        %
        for i=1:2
            text4 = fgetl(fid);
        end
        noddata = fscanf(fid,'%i %f %f %f %f %f %f',[7 npoin]);
        noddata = noddata';
        %
        nodid = squeeze(noddata(:,1));
        coord(nodid,1) = squeeze(noddata(:,2));
        coord(nodid,2) = squeeze(noddata(:,3));
        unknp(nodid,1) = squeeze(noddata(:,4));
        unknp(nodid,2) = squeeze(noddata(:,5));
        unknp(nodid,3) = squeeze(noddata(:,6));
        unknp(nodid,4) = squeeze(noddata(:,7));
        %
        for ip = 1:npoin
            if noddata(ip,1) ~= ip
                fprintf(' node no. %5d in data file is missing', ip);
                pause;
                break
            end
        end
        %
        for i=1:2
            text5 = fgetl(fid);
        end
        eledata = fscanf(fid,'%i %i %i %i',[4 nelem]);
        eledata = eledata';
        %
        eleid = squeeze(eledata(:,1));
        intmat = [];
        intmat(eleid,1) = squeeze(eledata(:,2));
        intmat(eleid,2) = squeeze(eledata(:,3));
        intmat(eleid,3) = squeeze(eledata(:,4));
        %
        for ie = 1:nelem
            if eledata(ie,1) ~= ie
                fprintf(' element no. %5d in data file is missing', ie);
                pause;
                break
            end
        end
        %
        for i=1:2
            text6 = fgetl(fid);
        end
        boudata = fscanf(fid,'%i %i %i %i %i',[4 nboun]);
        bsido = boudata';
        %


        %  Compute total numer of sides.
        %
        nside = (3.*nelem + nboun)/2.;
        %
        fprintf('\n');
        fprintf('The finite element model consists of: \n');
        fprintf('     number of nodes                   = %7d\n',npoin);
        fprintf('     number of elements                = %7d\n',nelem);
        fprintf('     number of boundary conditions     = %7d\n',nboun);
        fprintf('     number of iterations needed       = %7d\n',niter);
        fprintf('\n');
        %
        %  Save these nodal initial conditions:
        %
        unknp1 = unknp;
        %
        %  Obtain nodal conservation variables:
        %
        for j=2:4
            for i=1:npoin
                unknp(i,j) = unknp(i,1)*unknp(i,j);
            end 
        end
        %
        %  Compute element quantities:
        %
        unkno = zeros(nelem,4);
        factor = 1/3; ncount = 3;
        for ie=1:nelem
            for ia=1:4
                unkno(ie,ia) = 0.;
                for in=1:ncount
                    ip = intmat(ie,in);
                    unkno(ie,ia) = unkno(ie,ia) + factor*unknp(ip,ia);
                end
            end
        end
        %
        %  Identify the side: determine arrays iside(nside,4), jesid(nelem,4) and 
        %  rsido(nside,3):
        %
        rsido = zeros(nside,3); lhowm = zeros(npoin); jesid = zeros(nelem,3);
        %
        for ie=1:nelem
            ncount = 3;
            for in=1:ncount
                ip = intmat(ie,in);
                lhowm(ip) = lhowm(ip) + 1;
            end
        end
        %
        %  Fill in array lwher: location of each node inside icone:
        %
        lwher = zeros(npoin,1); 
        lwher(1) = 0;
        for ip=2:npoin
            lwher(ip) = lwher(ip-1) + lhowm(ip-1);
        end 
        %
        %  Fill in array icone: elements in each node:
        %
        lhowm = zeros(npoin,1); icone = zeros(4*nelem,1);
        for ie=1:nelem
            ncount = 3;
            for in=1:ncount
                ip = intmat(ie,in);
                lhowm(ip) = lhowm(ip) + 1;
                jloca = lwher(ip) + lhowm(ip);
                icone(jloca) = ie;
            end 
        end
        %
        %  Loop over the nodes:
        %
        iside = zeros(nside,4); 
        iloca = 0;
        for ip=1:npoin
            iloc1 = iloca;
            iele = lhowm(ip);
            if iele~=0
                iwher = lwher(ip);
        %
        %  Loop over elements surrounding the point ip:
        %    
                ip1 = ip;
                for iel=1:iele
                    ie = icone(iwher+iel);
        %  
        %  Find out position of ip in the conectivity matrix:
        %
                    ncount = 3;
                    for in=1:ncount
                        in1 = in;
                        ipt = intmat(ie,in);
                        if ipt==ip 
                            break
                        end
                    end  
                    j = 0;            
                    for jnod=1:ncount-2:ncount-1
                        j = j + 1;
                        in2 = in1 + jnod;
                        if in2>ncount
                            in2 = in2 - ncount;
                        end
                        ip2 = intmat(ie,in2);               
                        if ip2>=ip1
                            i = 0;
        %                    
        % Check the side: new or old
        %
                            if iloca~=iloc1
                                for is=iloc1+1:iloca
                                    jloca = is;
                                    if iside(is,2)==ip2
                                        iside(jloca,2+j) = ie;
                                        i = 1;
                                        break
                                    end
                                end
                            end
                            if i==0
                                iloca = iloca + 1;
                                iside(iloca,1) = ip1;
                                iside(iloca,2) = ip2;
                                iside(iloca,2+j) = ie;
                            end
                        end
                    end
        %
        %  End loop over elements surrounding point ip:
        %
                end
        %
                for is=iloc1+1:iloca
                    if iside(is,3)==0
                        iside(is,3) = iside(is,4);
                        iside(is,4) = 0;
                        iside(is,1) = iside(is,2);
                        iside(is,2) = ip1;
                    end
                end
        %
        %  End loop over points
        %
            end
        end
        %
        %  Now reset the boundary markers:
        %
        for is=1:nside
            if iside(is,4)==0
                il = iside(is,1);
                ir = iside(is,2);
                ie = iside(is,3);
                for ib=1:nboun
                    ibe = bsido(ib,3);
                    if ibe==ie
                        ilb = bsido(ib,1);
                        irb = bsido(ib,2);  
                        if ilb==il && irb==ir
                            iside(is,4) = -bsido(ib,4);
                            break
                        end
                    end
                end
            end
        end
        %
        %  Form the element/sides connectivity array:
        %
        for is=1:nside
            iel = iside(is,3);
            ier = iside(is,4);
            inode = iside(is,1);
            jnode = iside(is,2);
            ncount = 3;
            for in=1:ncount
                i1  = intmat(iel,in);
                in1 = in + 1;
                if in1>ncount
                    in1 = 1;
                end
                i2 = intmat(iel,in1);
                if (inode==i1) && (jnode==i2)
                    jesid(iel,in) = is;
                end
            end
            if ier>0
                ncount = 3;
                for in=1:ncount
                    i1 = intmat(ier,in);
                    in1 = in + 1;
                    if in1>ncount
                        in1 = 1;
                    end
                    i2 = intmat(ier,in1);
                    if (inode==i2) && (jnode==i1) 
                        jesid(ier,in) = is;
                    end
                end
            end
        end
        %
        %  Compute components of unit normal vector to the side and the lengths:
        %
        for is=1:nside 
            ipi = iside(is,1);
            ipj = iside(is,2);
            dx = coord(ipj,1) - coord(ipi,1);
            dy = coord(ipj,2) - coord(ipi,2);
            dl = sqrt(dx^2 +dy^2);
            rsido(is,1) =  dy/dl;
            rsido(is,2) = -dx/dl;
            rsido(is,3) =  dl;
        end
        %
        %  Compute element areas:
        %
        x = zeros(3,1); y = zeros(3,1); area = zeros(nelem,1);  
        for ie=1:nelem
            for ia=1:3
                n = intmat(ie,ia);
                x(ia) = coord(n,1); y(ia) = coord(n,2);
            end
            b1 = y(2) - y(3); b2 = y(3) - y(1); b3 = y(1) - y(2);
            c1 = x(3) - x(2); c2 = x(1) - x(3); c3 = x(2) - x(1);
            area(ie) = 0.5*( x(1)*b1 + x(2)*b2 + x(3)*b3 );
        end
        %
        %  Determine representative 'element lengths' for time step computation.
        %  There are 2 values for each side, eah pepresents the normal distance 
        %  from the element centroids(2 elements on both sides) to the side 
        %  considered.
        %
        tol = 1.0e-10;
        %
        %  Loop over number of elements:
        %
        slen = zeros(nside,2); delte = zeros(nelem,1); 
        for ie=1:nelem
            xc = 0.0; 
            yc = 0.0;
            for in=1:3
                ip = intmat(ie,in);
                xc = xc + coord(ip,1); 
                yc = yc + coord(ip,2);
            end
            xc = xc/3; 
            yc = yc/3;
            for in=1:3
                is = jesid(ie,in);
                il = iside(is,1);
                ir = iside(is,2);
                if il==intmat(ie,in) 
                    ies = 1;
                else
                    ies = 2;
                end
                xl = coord(il,1); 
                yl = coord(il,2);
                xr = coord(ir,1); 
                yr = coord(ir,2);
                if abs(yr-yl)<tol
                    dist = abs(yc-yl);
                else
                    if abs(xr-xl)>tol
                        rm = (yr-yl)/(xr-xl);
                        c = yl - rm*xl;
                        rm1 =  -1.0/rm;
        %
        % Compute intersection point:
        %
                        xp =(yc - rm1*xc - c)/(rm - rm1);
                        yp = rm1*(xp - xc) + yc;
                        dist = sqrt( (xp - xc)^2 + (yp - yc)^2 );
                    else
                        dist = abs(xc - xl);
                    end
                end
                slen(is,ies) = dist;
        %        
        %  Deal with boundary element lengths:
        %
                iel = iside(is,3);
                ier = iside(is,4);
                if iel <= 0
                    slen(is,1) = slen(is,2);
                else
                    if ier <= 0
                        slen(is,2) = slen(is,1);
                    end
                end
            end
        end
        %
        ishow = 1;
        fprintf('  Performing iterations for convergence  \n');
        fprintf('   Iter      Del rho    Del rho-u    Del rho-v    Del rho-e\n');
        %
        %  Enter iteration loop:
        %
        % Setup iteration for plot
        % Setup residual plotting animation
        h1 = animatedline;
        h1.Color = 'black';
        h1.DisplayName = 'Density';

        h2 = animatedline;
        h2.Color = 'blue';
        h2.DisplayName = 'X-Velocity';

        h3 = animatedline;
        h3.Color = 'green';
        h3.DisplayName = 'Y-Velocity';

        h4 = animatedline;
        h4.Color = 'red';
        h4.DisplayName = 'Energy';
        set(gcf,'currentchar',' ');
        set(gca, 'YScale', 'log');
        xp = 1:niter*10;
        sumsqp = zeros(niter*10,4);
        % Figure for loop stopping
        stop = figure(1);
        iter = 0;
        while 1
           iter = iter + 1;  %  ***  Start iteration loop  ***
            %
            %  Store element unknowns of previous iteration in unkn1(nelem,4)
            %
            unkn1 = unkno;
            gam1 = gamma - 1;
            %    
            %  Initialize element time steps:
            %
                for ie=1:nelem
                    delte(ie) = 1.0e+10;
                end
            %
            %  Initialze rhs0 vector:
            %
                rhs0 = zeros(nelem,4); rlam = zeros(4,1);
                diss = zeros(4,1); flux = zeros(4,1); sumsq = zeros(4,1); 
                r = zeros(4,4); avroe = zeros(4,4);
            %
            %  Loop over the sides:
            %
            for is=1:nside
            %
            %  Identify the left and right element numbers:
            %
            iel = iside(is,3);
            ier = iside(is,4);
                if ier==0 
                    fprintf('\n');
                    fprintf('\n         !!  Data Error  !!        ');
                    fprintf('\n *** Please Check Boundary Data ***');
                    fprintf('\n');
                    return
                end
            %
            %  Get components of normal vector for the side considered:
            %
                rnx = rsido(is,1);
                rny = rsido(is,2);
                rlen= rsido(is,3);
            %
            %  collect the "left" element values:
            %
                rhol = unkno(iel,1);
                uxl  = unkno(iel,2)/rhol;
                vyl  = unkno(iel,3)/rhol;
                tel  = unkno(iel,4)/rhol;
                presl= gam1*( unkno(iel,4) - 0.5*rhol*(uxl^2 + vyl^2) );
            %
            %  "left" normal and tangential velocities and total enthalpy:
            %
                unl =  uxl*rnx + vyl*rny;
                vtl = -uxl*rny + vyl*rnx;
                ul2 =  unl^2 + vtl^2;
                hl  =  gamma*tel - 0.5*gam1*ul2;
            %
            %  Is this side on the actual flow boundary ?
            %
            %  The right side is connected to actual element:
            %
                if ier>0
                    rhor = unkno(ier,1);
                    uxr  = unkno(ier,2)/rhor;
                    vyr  = unkno(ier,3)/rhor;
                    ter  = unkno(ier,4)/rhor;
                    presr= gam1*( unkno(ier,4) - 0.5*rhor*(uxr^2 + vyr^2) ); 
            %   
            %  Supersonic inflow:
            %
                elseif ier==-1
                    ii = iside(is,1);
                    jj = iside(is,2);
                    rhor = 0.5*( unknp(ii,1) + unknp(jj,1) );
                    uxr  = 0.5*( unknp(ii,2)/unknp(ii,1) + unknp(jj,2)/unknp(jj,1) );
                    vyr  = 0.5*( unknp(ii,3)/unknp(ii,1) + unknp(jj,3)/unknp(jj,1) );
                    ter  = 0.5*( unknp(ii,4)/unknp(ii,1) + unknp(jj,4)/unknp(jj,1) );
                    presr= gam1*( rhor*ter - 0.5*rhor*(uxr^2 + vyr^2) );
            % 
            %  Supersonic outflow:
            %
                elseif ier==-2
                    rhor = rhol;
                    uxr  = uxl;
                    vyr  = vyl;
                    ter  = tel;
                    presr= presl;
            %   
            %  Inviscid wall:
            %
                elseif ier==-3
                    rhor =  rhol;
                    uxr  = -rnx*unl - rny*vtl;
                    vyr  = -rny*unl + rnx*vtl;
                    presr=  presl;
                    ter  = (presr/(gam1*rhor)) + 0.5*(uxr^2 + vyr^2 );
                end
            %
            %  "right" normal and tangential velocities and total enthalpy:
            %
                unr =  uxr*rnx + vyr*rny;
                vtr = -uxr*rny + vyr*rnx;
                ur2 =  unr^2 + vtr^2;
                hr  =  gamma*ter - 0.5*gam1*ur2;
            %
            %  Compute interface values:
            %
                bi = sqrt(rhor/rhol);
                ai = 1/(1 + bi);
                ui = (bi*uxr + uxl)*ai;
                vi = (bi*vyr + vyl)*ai;
                hi = (bi*hr + hl)*ai;
                ci2= gam1*( hi - 0.5*(ui^2 + vi^2) );
                if ci2<0
                    fprintf('\n Negative sound speed between elements %5d %5d\n', ier, erel);
                    boom
                end
                ci   =  sqrt(ci2);
                ucap =  ui*rnx + vi*rny;
                vcap = -ui*rny + vi*rnx;
                cx   =  ci*rnx;
                cy   =  ci*rny;
                alp  =  0.5*(ui^2 + vi^2);
            %
            %  Compute the four absolute eigenvalues:
            %
                rlam(1) = abs(ucap);
                rlam(2) = abs(ucap);
                rlam(3) = abs(ucap + ci);
                rlam(4) = abs(ucap - ci);
            %
            %  Reset these eigenvalues so that the range is from zero to one:
            %
                eigmax = abs(ucap) + ci;
                for ir=1:4
                    rlam(ir) = rlam(ir)/eigmax;
                end
            %
            %  Set epslam = 0.01 (see reference [3] of chaper 10):
            %
                epslam =0.01;
                for ir=1:4
                    if rlam(ir)>=epslam
                        break
                    end
                    rlam(ir) = 0.5*(rlam(ir)^2/epslam + epslam);
                end
            %
            %  Reset back the correct (dimension) eigenvalues:
            %
                for ir=1:4
                    rlam(ir) = rlam(ir)*eigmax;
                end
            %
            %  Compute element time step associated with this side:
            %
                replen = slen(is,1) + slen(is,2);
                eigmax = abs(ucap) + ci;
            %
            %  Given Time step safety factor = 0.5:
            %

                dtl = csafe*replen/eigmax;
                delte(iel) = min(delte(iel),dtl);
                if ier>0
                    delte(ier) = min(delte(ier),dtl);
                end
            %
            %  Compute [r] matrix:
            %
                r(1,1) =  alp*gam1 - ci2;
                r(1,2) = -gam1*ui;
                r(1,3) = -gam1*vi;
                r(1,4) =  gam1;
                r(2,1) = -vcap;
                r(2,2) = -rny;
                r(2,3) =  rnx;
                r(2,4) =  0;
                r(3,1) =  alp*gam1 - ucap*ci;
                r(3,2) =  cx - gam1*ui;
                r(3,3) =  cy - gam1*vi;
                r(3,4) =  gam1;
                r(4,1) =  alp*gam1 + ucap*ci;
                r(4,2) = -cx - gam1*ui;
                r(4,3) = -cy - gam1*vi;
                r(4,4) =  gam1;
            %
            %  Compute [r] matrix inverse:
            %
                ri(1,1) = -1/ci2;
                ri(1,2) =  0;
                ri(1,3) =  0.5/ci2;
                ri(1,4) =  0.5/ci2;
                ri(2,1) = -ui/ci2;
                ri(2,2) = -rny;
                ri(2,3) = (ui + cx)/(2*ci2);
                ri(2,4) = (ui - cx)/(2*ci2);
                ri(3,1) = -vi/ci2;
                ri(3,2) =  rnx;
                ri(3,3) = (vi + cy)/(2*ci2);
                ri(3,4) = (vi - cy)/(2*ci2);
                ri(4,1) = -alp/ci2;
                ri(4,2) =  vcap;
                ri(4,3) = (alp + ucap*ci)/(2*ci2) + 1/(2*gam1);
                ri(4,4) = (alp - ucap*ci)/(2*ci2) + 1/(2*gam1);
            %
            %  Compute [as] = [ri] [eig] [r]:
            %
                for i=1:4
                    for j=1:4
                        r(i,j) = rlam(i)*r(i,j);
                    end
                end
            %
                for i=1:4
                    for j=1:4
                        avroe(i,j) = 0;
                            for l=1:4
                                avroe(i,j) = avroe(i,j) + ri(i,l)*r(l,j);
                            end
                    end
                end
            %
            %  Compute the difference of the conservation variables between 
            %  the right and left elements:
            %
                du(1) = rhol - rhor;
                du(2) = rhol*uxl - rhor*uxr;
                du(3) = rhol*vyl - rhor*vyr;
                du(4) = rhol*tel - rhor*ter;
            %
            %  Compute diss = [as] ul - ur:
            %
                for i=1:4
                    diss(i) = 0;
                        for j=1:4
                            diss(i) = diss(i) + avroe(i,j)*du(j);
                        end
                end
            %
            %  Compute sum of the left and the right fluxes:
            %
                fsum(1) = rhol*unl + rhor*unr;
                fsum(2) = rnx*(presl + presr) + rhol*uxl*unl + rhor*uxr*unr;
                fsum(3) = rny*(presl + presr) + rhol*vyl*unl + rhor*vyr*unr;
                fsum(4) = (rhol*tel + presl)*unl + (rhor*ter + presr)*unr;
            %
            %  The inviscid flux on rhs of the eq. is:
            %
                for i=1:4
                    flux(i) = 0.5*(fsum(i) + diss(i));
                end
            %
            %  Contribution of this flux to the "left" element:
            %
                for i=1:4
                    rhs0(iel,i) = rhs0(iel,i) - rlen*flux(i);
                end
            %
            %  Contribution of this flux to the "right" element:
            %
                if ier>=0
                    for i=1:4
                        rhs0(ier,i) = rhs0(ier,i) + rlen*flux(i);
                    end
                end
            % 
            %  End loop over all the sides:
            %
            end
            %
            %  Solve for nodal increments and update conservation variables:
            %
                for ie=1:nelem
                    for ia=1:4
                        unkno(ie,ia) = unkn1(ie,ia) + delte(ie)*rhs0(ie,ia)/area(ie);
                    end
                end
            %
            %  Chack for convergence:
            %
                for ia=1:4
                    sumsq(ia) = 0.;
                    for ie=1:nelem
                        diff = unkno(ie,ia) - unkn1(ie,ia);
                        sumsq(ia) = sumsq(ia) + diff^2;
                    end
                    sumsq(ia) = sqrt(sumsq(ia));
                end
                sumsqp(iter,:) = real(sumsq(:));
            %
                addpoints(h1,xp(iter),sumsqp(iter,1));
                addpoints(h2,xp(iter),sumsqp(iter,2));
                addpoints(h3,xp(iter),sumsqp(iter,3));
                addpoints(h4,xp(iter),sumsqp(iter,4));
                title("Press F to stop the iteration loop")
                legend
                drawnow
                pause(0);

                if iter==ishow
                    fprintf('%7d %10.5e %10.5e %10.5e %10.5e\n',iter, real(sumsq(1)), real(sumsq(2)), real(sumsq(3)), real(sumsq(4)));
                    if ishow~=1
                        ishow = ishow + nshow;
                    else
                        ishow = ishow + nshow - 1;
                    end
                end
                 if  get(gcf,'currentchar') == 'f'
                     break
                 end
                 if iter == niter

                    continueCalculation = false
                    while true

                        addniterBoolean = input("Do you want to continue calculating? (Y/N) :\n",'s');

                        % handle wrongly typed input
                        if addniterBoolean(1) ==  'y' || addniterBoolean(1) == 'Y'
                            addniter = input("How many more iteration?\n",'s');
                            addniter = str2double(addniter);
                            niter = niter + addniter
                            continueCalculation = true
                            break
                        elseif addniterBoolean(1) ==  'n' || addniterBoolean(1) == 'N'
                            break
                        else
                            continue
                        end
                    end

                    if (continueCalculation == false)
                        break
                    end
                    
                     
                 end
        end
        %
        %  ***  End iteration loop  ***
        %
        %  Print out solutions of density, u, v, velocity and total energy
        %
        unknp = zeros(npoin,4); idum1 = zeros(npoin);
        %
        for ie=1:nelem
            for in=1:3
                ip = intmat(ie,in);
                idum1(ip) = idum1(ip) + 1;
                for ia=1:4
                    unknp(ip,ia) = unknp(ip,ia) + unkno(ie,ia);
                end
            end
        end
        %
        for ia=1:4
            for ip=1:npoin
                unknp(ip,ia) = unknp(ip,ia)/idum1(ip);
            end
        end
        %
        %  Transform conservation variables back to primative variables:
        %
        for ip=1:npoin
            for ia=2:4
                unknp(ip,ia) = unknp(ip,ia)/unknp(ip,1);
            end
        end
        %
        %  Constraint some nodal quantities to inlet boundary conditions:
        %
        for ib=1:nboun
            ibc = bsido(ib,4);
            if ibc==1
                ii = bsido(ib,1);
                jj = bsido(ib,2);
                for ia=1:4
                    unknp(ii,ia) = unknp1(ii,ia);
                    unknp(jj,ia) = unknp1(jj,ia);
                end
            end
        end
        %
        unkno1 = zeros(npoin,1); unkno2 = zeros(npoin,1); unkno3 = zeros(npoin,1);
        unkno4 = zeros(npoin,1); unkno5 = zeros(npoin,1); 
        for ip=1:npoin
            unkno1(ip) = unknp(ip,1); 
            unkno2(ip) = unknp(ip,2); 
            unkno3(ip) = unknp(ip,3); 
            unkno4(ip) = unknp(ip,4); 
            unkno5(ip) = gam1*unkno1(ip)*( unkno4(ip) - 0.5*(unkno2(ip)^2 + unkno3(ip)^2) );
        end
        %
        %  Print out the nodal solutions into a file and on the screen:
        %
        filename = input('\n Enter file name for the solutions: ', 's');
        if filename ==  ' '
            filename = 'solution';
        end
        fid = fopen(filename, 'w');
        fprintf(fid,' Nodal density, velocity and total energy solutions [%5d]:\n', npoin);
        fprintf(    ' Nodal density, velocity and total energy solutions [%5d]:\n', npoin);
        fprintf(fid,'\n   Node        density        u-velocity        v-velocity         total energy\n');
        fprintf(    '\n   Node        density        u-velocity        v-velocity         total energy\n');

        for ip = 1:npoin
            fprintf(fid,' %6d  %16.6e  %16.6e  %16.6e  %16.6e\n', ip, unkno1(ip), unkno2(ip), unkno3(ip), unkno4(ip));
            fprintf(    ' %6d  %16.6e  %16.6e  %16.6e  %16.6e\n', ip, unkno1(ip), unkno2(ip), unkno3(ip), unkno4(ip));
        end
        %
        %  Prepare data for plotting on screen:
        %
        xmin = 999999; ymin = 999999; xmax = -999999; ymax = -999999;
        x = zeros(npoin,1); y = zeros(npoin,1);
        for i=1:npoin
            x(i)=coord(i,1);
            y(i)=coord(i,2);
        end
        for i=1:npoin
            if x(i) < xmin
                xmin = x(i);
            end
            if y(i) < ymin
                ymin = y(i);
            end
            if x(i) > xmax
                xmax = x(i);
            end
            if y(i) > ymax
                ymax = y(i);
            end
        end

    catch
        safe = input('Error detected in hiflow, (r)educe relaxation factor or (s)top the calculation\n', 's');
        if safe == 'r' || safe == 'R'
            if csafe < 0.1
                csafe = csafe*0.9
            else
                csafe = csafe - 0.05
            end
            if csafe < 0.0001
                fprintf('Try remeshing using new parameter \n')
                return
            end
            pause(2);
        elseif safe == 's' || safe == 'S'
            fprintf('Stopping the calculation, please use new input data')
            return
        end
        continue
    end
    
    %
    %  Plot the mesh with nodes and element numbers on the screen if needed:
    %
    clf

    while (true)
        a1 = input('Would you like to see the mesh plot? y/n: ', 's');
        
        % handle wrongly typed input
        if(a1(1) =='y')||(a1(1) =='Y');
            %
            %  Plot elements with node and element numbers:
            %
            TRI = [intmat(:,1) intmat(:,2) intmat(:,3)];
            X = coord(:,1); Y = coord(:,2); Z = 0.*X;
            figure(2+(cv-1)*6)
            trimesh(TRI,X,Y,Z,'edgecolor','k')
            view(2),axis([xmin xmax ymin ymax]),axis equal, grid off, hold on
            %    
            %  Place element numbers:
            %
            xe = zeros(nelem,1); ye = zeros(nelem,1);
            break
        elseif a1(1) == 'n' || a1(1) == 'N'        
            break
        else
            continue
        end
    end


    %  Display the density fringe plot on the screen if needed:
    %

    while (true)
        a2 = input('Would you like to see the density fringe plot? y/n: ', 's');
        
        % handle wrongly typed input
        if(a2(1) =='y')||(a2(1) =='Y');
            
    
            figure(3)
            trisurf(intmat,x,y,0*x,unkno1,'edgecolor','none','facecolor','interp');
            view(2),axis([xmin xmax ymin ymax]),axis equal,colorbar, grid off
        
            break
        elseif a2(1) == 'n' || a2(1) == 'N'        
            break
        else
            continue;
        end
    end
    
    
    %

    while (true)
        a3 = input('Would you like to see the velocity vectors plotted on mesh? y/n: ', 's');
        
        % handle wrongly typed input
        
        if(a3(1) =='y')||(a3(1) =='Y')
            %
            %  Perform vector plot on screen if needed:
            %  Plot elements:
            %
            TRI = [intmat(:,1) intmat(:,2) intmat(:,3)];
            X = coord(:,1); Y = coord(:,2); Z = 0.*X;
            figure(4)
            trimesh(TRI,X,Y,Z,'edgecolor','g')
            view(2),axis([xmin xmax ymin ymax]),axis equal, grid off, hold on
            %    
            xx = zeros(npoin,1); yy = zeros(npoin,1); uu = zeros(npoin,1); vv = zeros(npoin,1);
            for ip = 1:npoin
                xx(ip) = coord(ip,1); yy(ip) = coord(ip,2);
                uu(ip) = unkno2(ip); vv(ip) = unkno3(ip);
            end
            quiver(xx, yy, uu, vv,'r')
            break
        elseif a3(1) == 'n' || a3(1) == 'N'        
            break
        else
            continue
        end
    end

    %

    while (true)
        a4 = input('Would you like to see the pressure fringe plot? y/n: ', 's');
        
        % handle wrongly typed input
        if(a4(1) =='y')||(a4(1) =='Y');
            %
            %  Perform pressure fringe plot on screen if needed:
            %
            figure(5)
            trisurf(intmat,x,y,0*x,unkno5,'edgecolor','none','facecolor','interp');
            view(2),axis([xmin xmax ymin ymax]),axis equal, colorbar, grid off
            break
        elseif a4(1) == 'n' || a4(1) == 'N'        
            break
        else
            continue
        end
    end
    
    

    while (true)
        a5 = input('Would you like to see the residual plot? y/n: ', 's');
        
        % handle wrongly typed input
        
        if(a5(1) =='y')||(a5(1) =='Y');
            %
            %  Perform pressure fringe plot on screen if needed:
            %
            figure(6)
    
            xpf = xp(1:iter);
            % Repopulate residual plot
            sumsqpf = zeros(iter,4);
            plot(xpf, sumsqp(1:iter,1),'Color','black','DisplayName','Density')
            hold on
            plot(xpf, sumsqp(1:iter,2),'Color','blue','DisplayName','X-Velocity')
            hold on
            plot(xpf, sumsqp(1:iter,3),'Color','green','DisplayName','Y-Velocity')
            hold on
            plot(xpf, sumsqp(1:iter,4),'Color','red','DisplayName','Energy')
            set(gca, 'YScale', 'log');
            legend
            break
        elseif a5(1) == 'n' || a5(1) == 'N'        
            break
        else
            continue
        end
    end
    
    
    while (true)
        re = input('Would you like to apply adaptive remeshing? y/n:','s');
        
        % handle wrongly typed input
        if(re(1) =='y')||(re(1) =='Y')
            cv_s = int2str(cv);
            n1body = 1;
            % count for wall boundary condition for remeshing
            boudata = boudata';
            n2body = 0;
            fixm = zeros(nboun,4);
            fixh = zeros(nboun,2);
            intmat_bound = zeros(nboun,3);
            coor_bound = [nboun,6];
            coor_bound_local =[];
            for ibody = 1:nboun
    
                if boudata(ibody,4) == 3
                    n2body = n2body + 1;
                end
    
                fixm(ibody,1) = ibody;
                fixm(ibody,2) = 2;
                fixm(ibody,3) = 0;
                fixm(ibody,4) = boudata(ibody,4);
                fixh(ibody,1) = ibody;
                fixh(ibody,2) = ibody+1;
    
                if ibody == nboun
                    fixh(ibody,2) = 1;
                end
            end
    
    
            %% write .RE file for cv == 1
            cv = cv + 1;
            cv_s = int2str(cv);
            filetypeRE = fopen(strcat(filenameGEO,'.RE', cv_s),'w');
            fprintf(filetypeRE,'  NPOIN    NELEM       N1       N2 \n');
            fprintf(filetypeRE,'%d %d %d %d \n', npoin, nelem,n1body,n2body );
    
            for i = 1:npoin
    
                fprintf(filetypeRE,'%d %f %f %f %f %f %f \n', i, unkno1(i), unkno2(i), unkno3(i), unkno4(i), coord(i,1),coord(i,2));
    
            end
    
            for i = 1:nelem
                fprintf(filetypeRE,'%d %f %f %f \n', i, intmat(i,1), intmat(i,2), intmat(i,3));
            end
    
        % increase cv
            while 1
                % special case for only .dat
                % seperately generate .G1 file
                if (inp(1) == 'd')
                    % .G generator 
                    % .dat file generator
                    delkmma = 0;
                    delkmmi = 0;
                    delsca = 0;
                    cv;
                    cvnext_s = int2str(cv);
                    filetypeG = fopen(append(filenameDAT,'.G',cvnext_s),'w');
                    fprintf(filetypeG,'1 \n');
                    fprintf(filetypeG,'%f %f %f \n',delkmma*2,delkmmi/2,delsca/2);
                    fprintf(filetypeG,'%d \n',cv);
                    fprintf(filetypeG,'%d \n %d \n %d \n %d \n %d \n %d \n %d \n %d \n %d \n %d \n %d \n %d \n',1, 2, 5, 4, 3, 2, 3, 3, 4, 6, 5, 0);
                end

                % Call Adapt
                try
                    Adapt(filenameGEO,cv);
                catch MExc
                    disp(MExc)
                    fprintf('Please change the input value \n');
                    continue
                end
                
                resize = input('Do you want to resize the element again? (Y/N)\n','s');
                if resize(1) == 'y' || resize == 'Y'
                else
                    break
                end
            end
            break
        elseif re(1) == 'n' || re(1) == 'N'        
            return
        else
            continue
        end
    end
end
