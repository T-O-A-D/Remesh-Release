function [unkng, coorg, ieleg, delta, temp,geomg,fid_7,deri2,delkmma,delkmmi,delsca] = INPUT(neleg,npoig,ieleg,coorg,delta,unkng,unit_file,unit_7,io)

    disp(unit_file)
    delkmma = [];
    delkmmi = [];
    delsca = [];
    % extract each lines from file
    fid = fopen(unit_file, 'r');

    % ignore first 2 lines
    fgets(fid);
    fgets(fid);

    % get unkng and coorg
    temp = fscanf(fid,'%i %f %f %f %f %f %f',[7 npoig]);
    %temp = temp';
    unkng = temp(2:5,:);
    coorg = temp(6:7,:);

    % Single everything out
    unkng = single(unkng);
    
    coorg = single(coorg);
    
    % get ieleg
    temp = fscanf(fid,'%i %f %f %f',[4 neleg]);
    %temp = temp';
    ieleg = temp(2:4,:);

    % Single ieleg
    ieleg = single(ieleg);

    % GET THE VALUE OF THE REFINEMENT INDICATORS 
    %io = 0;
    %io = 1;
    if (io == '1')
        fid = fopen(unit_file, 'r');
        fgets(fid);
        npoig = fscanf(fid,'%f',1);
        neleg = fscanf(fid,'%f',1);
        n1 = fscanf(fid,'%f',1);
        n2 = fscanf(fid,'%f',1);
        coord = zeros(npoig,2);
        unkng = zeros(npoig,4);
        temp = fscanf(fid,'%i %f %f %f %f %f %f',[7 npoig]);
        temp = temp';
        unkng(:,1) = temp(:,2);
        unkng(:,2) = temp(:,3);
        unkng(:,3) = temp(:,4);
        unkng(:,4) = temp(:,5);
        coord(:,1) = temp(:,6);
        coord(:,2) = temp(:,7);
        % get ieleg
        temp = fscanf(fid,'%i %f %f %f',[4 neleg]);
        temp = temp';
        intmat = zeros(neleg,3);
        intmat(:,1) = temp(:,2);
        intmat(:,2) = temp(:,3);
        intmat(:,3) = temp(:,4);
        geomg = single(zeros(7,10000));
        nstop = single(0);
        intmat = single(intmat);
        coord = single(coord);
        unkng = single(unkng);
        [delta, unkng, geomg,fid_7,deri2,delkmma,delkmmi,delsca] = ref2(npoig,neleg,coord,geomg,nstop,intmat,unit_7,unkng);
    else
        %for ip = 1:npoig-1
        unkng(1,:) = input('Input initial uniform mesh size \n');
        for ip = 1:npoig
            delta(1,ip) = unkng(1,ip);
            delta(2,ip) = unkng(2,ip);
            delta(3,ip) = unkng(3,ip);
            delta(4,ip) = unkng(4,ip);
        end
        % READ THE VALUE OF THE UNKNOWNS FOR INTEPOLATION

        % extract each lines from file

        fid = fopen(unit_file, 'r');

        % ignore first all previous lines
        fgets(fid);
        fgets(fid);
        for ip = 1:npoig+neleg
            fgets(fid);
        end
        fgets(fid);

        %start scanning after INITIAL CONDITIONS
        temp = fscanf(fid,'%i %f %f %f %f',[5 npoig]);
        temp = temp';
        % Single temp
        temp = single(temp);
        unkng = temp(:,2:5);
        geomg = single([]);
        fid_7 = 0;
        deri2 = single([]);
    end
end