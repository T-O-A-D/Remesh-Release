function [delta, unkng, geomg, fid_7,deri2,delkmma,delkmmi,delsca] = ref2(npoig,neleg,coord,geomg,nstop,intmat,unit_7,unkng)

    % OBTAIN GLOBALLY THE GEOMETRICAL VALUES NEEDED 
    namat = 4;
    nnode = 3;
    ngeom = 7;
    coord = coord';
    unkng = unkng';
    intmat = intmat';
    % test section
    % call getgeo
    geomg = getgeo(neleg,npoig,nnode,ngeom,intmat,coord,geomg,nstop);

    % initialize
    for ipoin = 1:npoig
        delta(1,ipoin) = 1.e+6;
        delta(2,ipoin) = 1.e+6;
    end

    % INQUIRE FOR MAXIMUM AND MINIMUM SPACING
    % get DELKMMA,DELKMMI,DELSCA
    % extract each lines from file

    fid_7 = fopen(unit_7, 'r');
    fgetl(fid_7);
    temp = fscanf(fid_7,'%f %f %f',[3 1]);
    temp = temp';

    delkmma = input('Input maximum element size \n');
    delkmmi = input('Input minimum element size\n');
    delsca = input('Input checking element size\n');
    % ----------

    % get STR
    fprintf(' ENTER MAXIMUM STRETCHING \n')
    temp = fscanf(fid_7,'%f',[1 1]);
    temp = temp';
    str = temp(1);
    str = input('Input cell aspect ratio \n');
    % to single
    str = single(str);
    % ----------

    % INQUIRE FOR REFINEMENT CRITERIA
    disp(' WHAT DO YOU WANT TO USE AS KEY VARIABLE ?')
    disp(' 1 - DENSITY')
    disp(' 2 - VELOCITY (MODULUS)')
    disp(' 3 - DENSITY + VELOCITY')
    disp(' 4 - ENTROPY')
    disp(' 5 - DENSITY + MACH NO')
    disp(' 6 - MACH NUMBER')
    % get ICRIT
    temp = fscanf(fid_7,'%f',[1 1]);
    temp = temp';
    icrit = temp(1);
    icrit = input('Please input the key variable for adaptive remeshing :');
    disp(icrit)

    % to single
    icrit = single(icrit);
    % ----------

    iback=0;
    ivari=1;

    % to single
    iback = single(iback);
    ivari = single(ivari);
    % ----------

    if (icrit == 2 || icrit == 3)
        ivari=2;
    end
    if (icrit == 4)
        ivari=3;
    end
    if (icrit == 5 || icrit == 6)
        ivari=4;
    end
    if (icrit == 3 || icrit == 5)
        iback=1;
    end
    pass = true; % jump gate label 1210
    % OBTAIN THE GRADIENT AT THE NODES
    derip = single([]);
    deri2 = single([]);
    deria = single([]);
    mmatg = single(zeros(1,npoig));
    while true % loop 1209
        if pass == true
            ;
        else
            if (iback == 1)
                ivari=1;
                iback=0;
            end
        end
        % jump gate exit label 1210
         % call getvar
        [deri2 nx ny fid_7] = getvar(nnode,namat,ngeom,2,neleg,npoig,coord,intmat,derip,deri2,deria,geomg,mmatg,unkng,ivari,fid_7);
        % FIND OUT PRINCIPAL DIRECTIONS AND VARIATIONS
        for ipoin = 1:npoig % loop lebel 1601

            amean = 0.5*(deri2(1,ipoin)+deri2(4,ipoin));
            adevi = 0.5*(deri2(1,ipoin)-deri2(4,ipoin));
            tauxy = deri2(2,ipoin);
            si1   = amean+sqrt(adevi^2+tauxy^2);
            si2   = amean-sqrt(adevi^2+tauxy^2);

            if (abs(si1) + abs(si2) < 1e-08)

                help(1,ipoin) = 1.e-09;
                help(2,ipoin) = 1.e-09;
                help(3,ipoin) = 1.0;
                help(4,ipoin) = 0.0;

            else
                if (amean >= 0.0)

                    help(1,ipoin) = abs(si1);
                    help(2,ipoin) = abs(si2);
                    alpha         = atan2(tauxy,adevi)/2.0;
                    help(3,ipoin) = -sin(alpha);
                    help(4,ipoin) = cos(alpha);

                else

                    help(1,ipoin) = abs(si2);
                    help(2,ipoin) = abs(si1);
                    alpha         = atan2(tauxy,adevi)/2.0;
                    help(3,ipoin) = cos(alpha);
                    help(4,ipoin) = sin(alpha);

                end
            end
        end

        % SEARCH FOR THE MAXIMUM SPACING AND STRETCHING
        amax=0.0;
        for ipoin = 1:npoig % loop label 1201
            if (help(1,ipoin) > amax)
                amax = help(1,ipoin);
            end
        end

        amaxsq=sqrt(amax);
        for ipoin = 1:npoig % loop label 1202

            delkm=delkmmi*amaxsq/sqrt(amax1([1.e-9,help(1,ipoin)]));
            delt2=delkmmi*amaxsq/sqrt(amax1([1.e-9,help(2,ipoin)]));

            if (delkm > delkmma)
                delkm = delkmma;
            end

            if (delt2 > delkmma)
                delt2 = delkmma;
            end

            if (delkm < delsca)
                delkm = delsca;
            end

            help(1,ipoin) = delkm;
            help(2,ipoin) = delt2;
        end

        for ipoin = 1:npoig % loop label 1302

            delkm = help(1,ipoin);
            delta(2,ipoin) = amin1([delta(2,ipoin),help(2,ipoin)]);

            if (delkm >= delta(1,ipoin))
                continue
            else
                delta(1,ipoin) = delkm;
                delta(3,ipoin) = help(3,ipoin);
                delta(4,ipoin) = help(4,ipoin);
            end

        end

        if (iback == 1) % jump gate label 1209
            pass = false;
            continue
        else
            break
        end
    end

    for ipoin = 1:npoig % loop label 1402
        delkm = delta(1,ipoin);
        delt2 = delta(2,ipoin);
        stre  = delt2/delkm;

        if (stre > str)
            stre = str;
        end
        if (stre < 1.0)
            stre = 1.0;
        end

        delta(2,ipoin) = stre;
    end

    % SOME EXTRA REFINEMENT FOR STAGNATION POINTS 
    disp(" DO YOU WANT EXTRA REFINEMENT FOR STAGNATION POINTS ?")

    % get ISTAG
    temp = fscanf(fid_7,'%f',[1 1]);
    temp = temp';
    istag = temp(1);

    % to single
    istag = single(istag);
    disp(istag)
    if (istag == 0)
        return
    end
    disp('ENTER DELKM, THRESHOLD MACH NUMBER AND GAMMA')
    % get delst,amini,gamma
    temp = fscanf(fid_7,'%f %f %f',[3 1]);
    temp = temp';
    delst = temp(1);
    amini = temp(2);
    gamma = temp(3);
    disp(temp)
    % to single
    delst = single(delst);
    amini = single(amini);
    gamma = single(gamma);
    % ------------

    for ipoin = 1:npoig % loop label 3430

        ro    = unkng(1,ipoin);
        u1    = unkng(2,ipoin);
        u2    = unkng(3,ipoin);
        en    = unkng(4,ipoin);
        vsq   = u1*u1 + u2*u2;
        pre   = (gamma-1.)*ro*(en-0.5*vsq);
        ropoi = sqrt((vsq*ro)/(gamma*pre));

        if (ropoi > amini)
            continue
        end

        delkm = delta(1,ipoin);
        delt2 = delkm*delta(2,ipoin);

        if (delkm > delst)
            delkm = delst;
        end

        delta(1,ipoin) = delkm;
    end

return
end