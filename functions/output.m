function output(npoig,neleg,coorg,intmag,intmeg,unkng,npoin,nelem,coord,intmat,rcond,nboun,lpoin,iside,nside,toler,ilast,filename,arZ,cv,delkmma,delkmmi,delsca)
    % this function is an output for input routine
    % read the unknowns for interpolation
    zero = 0;
    iz = 0;
    fileversion = cv;
    filnam = filename;
    unit_7  = append(filnam,'.G',   cv);
    unit_8  = append(filnam,'.KK', cv);
    unit_9  = append(filnam,'.RE',  cv);
    unit_10 = append(filnam,'.IN',  cv);
    unit_17 = append(filnam,'.FIX');
    unit_19 = append(filnam,'.SHO');
    unit_13 = append(filnam,'.D',  cv);
 

    % filetypeD = fopen(unit_13,'w');
    % for i = 1:npoin % execution loop label 100
    %     fprintf(filetypeD,'GRID           %5d       %d %7.4f %7.4f %7.4f       %d \n',i,iz,coord(1,i),coord(2,i),zero,iz);
    % end
    n = 1;
    l = 0;
    % for i = 1:nelem % execution loop label 200
    %     fprintf(filetypeD,'CTRIA3       %5d       %1d   %5d   %5d   %5d \n',i,n,intmat(1,i), intmat(2,i), intmat(3,i));
    % end 
    % filetypeKK = fopen(unit_8,'w');
    % % fprintf(filetypeKK,'NELEM NPOIN NBOUN\n');
    % % fprintf(filetypeKK,'%i %i %i\n',[nelem npoin nboun]);
    % 
    % filetypeIN = fopen(unit_10,'w');
    % fprintf(filetypeIN,'NELEM NPOIN NBOUN\n');
    % fprintf(filetypeIN,'%i %i %i\n',[nelem npoin nboun]);
    filetypeFIX = fopen(unit_17,'r');
    scanFIX = 0;
    while scanFIX == 0
        text = fgetl(filetypeFIX);
        if isequal(text,'INPUT GAMMA,EPSLAM,NTIME')
            scanFIX = 1;
        else
            scanFIX = 0;
        end
    end
    gamma = fscanf(filetypeFIX,'%f',1);
    epslam = fscanf(filetypeFIX,'%f',1);
    ntime = fscanf(filetypeFIX,'%f',1);
    temp = fgetl(filetypeFIX);
    temp = fgetl(filetypeFIX);
    p00 = fscanf(filetypeFIX,'%f',1);
    u00 = fscanf(filetypeFIX,'%f',1);
    v00 = fscanf(filetypeFIX,'%f',1);
    e00 = fscanf(filetypeFIX,'%f',1);
    % fprintf(filetypeIN,'INPUT GAMMA,EPSLAM,NTIME\n');
    % fprintf(filetypeIN,'%f %f %f\n',[gamma epslam ntime]);
    % fprintf(filetypeIN,'ELEMENT NODAL CONNECTIONS [ %i]:\n',nelem);
    % for i = 1:nelem % execution loop label 500
    %     fprintf(filetypeIN,'%i %i %i %i\n',[i intmat(1,i) intmat(2,i) intmat(3,i)]);
    % end
    % fprintf(filetypeIN,'NODAL COORDINATES [ %i]:\n',npoin);
    % for i = 1:npoin % execution loop label 510
    %     fprintf(filetypeIN,'%f %f %f\n',[i coord(1,i) coord(2,i) ]);
    % end
    % fprintf(filetypeIN,'NODAL INITIAL CONDITIONS [ %i]:\n',npoin);
    % fprintf(filetypeKK,'LAST VALUES FOR THE UNKNOWNS :\n');
    i1 = single([]);
    i2 = single([]);
    i3 = single([]);
    ienr = single([]);
    ilast = single(1);
    uv = u00;
    vv = v00;
    ro = p00;
    en = e00;
    %unkng = unkng';
    if cv == '1'
        unkng = unkng';
    end
    for i = 1:npoin % execution loop label 1002
        % interpolate the unknown
        inorm = 0;
        [iner,arZ, i1, i2, i3] = findel(npoig,neleg,coorg,intmag,intmeg,coord(1,i),coord(2,i),ilast,arZ,i1,i2,i3,ienr);
        inorm = 0;
        [dis, alp, anx, any] = getvalue(npoig,i1,i2,i3,arZ,unkng,p00,u00,v00,e00,inorm);
        if i > nboun % jump gate label 1202
        else
            if npoin > 0 % jump gate label 1202
            else
                tanx = rcond(2,i);
                tany = -rcond(1,i);
                veloc = tanx*uv + uany*vv;
                uv = veloc*tanx;
                vv = veloc*tany;
            end
        end
        % jump gate exit label 1202
    %     fprintf(filetypeKK,"%f %f %f %f %f %f %f \n",[i ro uv vv en coord(1,i) coord(2,i)]);
    %     fprintf(filetypeIN,'%f %f %f %f %f \n',[i p00 u00 v00 e00]);
    end

    % fprintf(filetypeKK,'INTERDEPENDENCY MATRIX INTMAT :\n');
    % for i = 1:nelem % execution loop label 1000
    %     fprintf(filetypeKK,'%i %i %i %i %i\n',[i intmat(1,i) intmat(2,i) intmat(3,i) 0]);
    % end
    % 
    % fprintf(filetypeKK,'BOUNDARY MARKERS AND NORMALS :\n');
    % 
    for i = 1:nboun % execution loop label 1003
        a = abs(rcond(1,i)) + abs(rcond(2,i));
        if a > 1*10^-5
    %         fprintf(filetypeKK,'%i %f %f %f\n',[i lpoin(i) rcond(1,i) rcond(2,i)]);
        elseif a < 1*10^-5
    %         fprintf(filetypeKK,'%i\n',[i]);
            lpoin(i) = 100;
        end
    end
    % 
    % fprintf(filetypeKK,'OUTPUT BSIDO \n');
    % fprintf(filetypeIN,'BOUNDARY CONDITIONS [ %i]: \n',[nboun]);
    % 
    for is = 1:nside % execution loop label 2500
        if iside(4,is) ~= 0 % jump gate label 2500
            continue
        else
            i1 = iside(1,is);
            i2 = iside(2,is);
            ie = iside(3,is);
            ib = lpoin(i1);
            if ib == 100
                ib = lpoin(i2);
            end
    %         fprintf(filetypeKK,' %i %i %i %i \n',[i1 i2 ie ib]);
    %         fprintf(filetypeIN,' %i %i %i %i\n',[i1 i2 ie ib]);
        end
    end
    % 


    % .dat file generator
    filenameDAT = append(filnam,'.dat',  cv);
    filetypeDAT = fopen(filenameDAT,'w');
    fprintf(filetypeDAT,'2 \n');
    fprintf(filetypeDAT,'FINITE ELEMENT MODEL FOR SHOCKWAVE \n');
    fprintf(filetypeDAT,'MODEL WITH %i TRIANGULAR ELEMENTS AND %i NODES \n', nelem, npoin);

    % write NELEM     NPOIN    NBOUN     NITER     NSHOW	to .dat
    fprintf(filetypeDAT,'NELEM     NPOIN    NBOUN     NITER     NSHOW	\n');
    fprintf(filetypeDAT,'%d %d %d %d %d \n', nelem, npoin, nboun, ntime, 10);

    % write NELEM     NPOIN    NBOUN     NITER     NSHOW	to .dat
    fprintf(filetypeDAT,' GAMMA	\n');
    fprintf(filetypeDAT,'%d \n', gamma);


    % read til initial condition
    fclose(filetypeFIX);
    filetypeFIX = fopen(strcat(filename,'.FIX'),'r');
    scanFIX = 0;
    while scanFIX == 0
        text = fgetl(filetypeFIX);
        if isequal(text,'INITIAL CONDITIONS')
            scanFIX = 1;
        else
            scanFIX = 0;
        end
    end

    % get initial condition
    init_1    = fscanf(filetypeFIX,'%f',1);
    init_2    = fscanf(filetypeFIX,'%f',1);
    init_3    = fscanf(filetypeFIX,'%f',1);
    init_4    = fscanf(filetypeFIX,'%f',1);

    % write coord and initial conditions to .dat
    fprintf(filetypeDAT,'NODAL COORDINATES AND INITIAL CONDITIONS [  ]:\n');
    for i = 1:npoin % execution loop label 100
        fprintf(filetypeDAT,'%d %f %f %f %f %f %f \n',i, coord(1,i),coord(2,i),init_1,init_2,init_3,init_4);
    end

    % write ELEMENT NODAL CONNECTIONS to .dat
    fprintf(filetypeDAT,' ELEMENT NODAL CONNECTIONS [   ]:\n');
    for i = 1:nelem % execution loop label 200
        fprintf(filetypeDAT,'%d %d %d  %d\n',i,intmat(1,i), intmat(2,i), intmat(3,i));
    end 

    % write ELEMENT NODAL CONNECTIONS FOR BOUNDARY CONDITIONS to .dat
    fprintf(filetypeDAT,' ELEMENT NODAL CONNECTIONS FOR BOUNDARY CONDITIONS [  ]:	\n');
    for is = 1:nside % execution loop label 2500
        if iside(4,is) ~= 0 % jump gate label 2500
            continue
        else
            i1 = iside(1,is);
            i2 = iside(2,is);
            ie = iside(3,is);
            ib = lpoin(i1);
            if ib == 100
                ib = lpoin(i2);
            end
            fprintf(filetypeDAT,' %i %i %i %i \n',[i1 i2 ie ib]);
        end
    end

    % .G generator 
    % .dat file generator

    cvnext = str2num(cv) + 1;
    cvnext_s = int2str(cvnext);
    filetypeG = fopen(append(filename,'.G',cvnext_s),'w');
    fprintf(filetypeG,'1 \n');
    fprintf(filetypeG,'%f %f %f \n',delkmma,delkmmi,delsca);
    fprintf(filetypeG,'%d \n',cvnext);
    fprintf(filetypeG,'%d \n %d \n %d \n %d \n %d \n %d \n %d \n %d \n %d \n %d \n %d \n %d \n',1, 2, 5, 4, 3, 2, 3, 3, 4, 6, 5, 0);
        
    return

end





