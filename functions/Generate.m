function Generate(filnamE)
%% PROGRAM REMESH
%  FINITE ELEMENT COMPUTER PROGRAM FOR GENERATING INITIAL AND ADAPTIVE MESH  
%    THE VALUES DECLARED IN THE PARAMETER STATEMENT BELOW SHOULD BE 
%    SET ACCORDING TO THE SIZE OF THE PROBLEMS :
%       MXPOI = MAXIMUM NUMBER OF NODES IN THE MODEL.
%       MXELE = MAXIMUM NUMBER OF ELEMENTS IN THE MODEL.
%       MXBOU = MAXIMUM NUMBER OF BOUNDARIES IN THE MODEL.

%% WRITE OUTP4UT FILE FOR ANSYS

filnam = filnamE;
% filnam = input(' PLEASE TYPE INPUT FILE NAME : ','s');
cv = "1";
%cv = '2';
% disp(filnam)
% disp(cv)


%% unit files
% UNIT  7 ==> INPUT DATA FROM SCREEN
% UNIT  8 ==> OUTPUT FILE 
% UNIT  9 ==> REMESH DATA FILE
% UNIT 10 ==> OUTPUT FILE FOR FINITE
% UNIT 17 ==> FIXED INPUT FILE
% UNIT 13 ==> INPUT NASTRAN FILE 
unit_7  = append(filnam,'.G',   cv);
unit_8  = append(filnam,'.KK', cv);
unit_9  = append(filnam,'.RE',  cv);
unit_10 = append(filnam,'.IN',  cv);
unit_17 = append(filnam,'.FIX');
unit_19 = append(filnam,'.SHO');
unit_13 = append(filnam,'.D',  cv);
disp('--------------------------------------------------')
%% define
% READ THE BACGROUND GRID DATA
lcoor = zeros(50000,1);
ieleg = zeros(3,10000);
nsidg = 0;
coorg = [];
%delta = zeros(4,5000);
delta = [];
unkng = [];
npfront = [];
nqfront = [];
nregi = [];
iside = [];
lwher = [];
lhowm = [];
icone = [];
intmeg = [];

igraph = 0;
ilast = 1;

%% main
ia = rfilliv(lcoor,5000,0);
disp(' *** READING DATA')

% CONTROL
unit_9_file = textread(unit_9,'%s', 'delimiter','\n','whitespace','');
texts = char(unit_9_file(1));
texts = allwords(texts);

%CONTROL
unit_9_secondline = char(unit_9_file(2));
unit_9_secondline = allwords(unit_9_secondline);
unit_9_secondline = str2double(unit_9_secondline);
npoig  = unit_9_secondline(1);
neleg  = unit_9_secondline(2);
n1body = unit_9_secondline(3);
n2body = unit_9_secondline(4);

%---------- to single ------------
npoig  = single(npoig );
neleg  = single(neleg );
n1body = single(n1body);
n2body = single(n2body);
delkmma = [];
delkmmi = [];
delsca = [];
% --------------------------------
io = '0';
[unkng, coorg, ieleg, delta, temp,geomg,fid_7,deri2,delkmma,delkmmi,delsca] = INPUT(neleg,npoig,ieleg,coorg,delta,unkng,unit_9,unit_7,io);
% ADD TRIANGLES TO MAKE THE DOMAIN COVERED BY THE BACKGROUND GRID CONVEX REGION
figure(1);
z = coorg(1,:).*0;
trimesh(ieleg',coorg(1,:)',coorg(2,:)',z,'EdgeColor','black');
view(2),axis equal, grid off, hold on
disp(' *** FILLING IN THE GAPS')

[npfront,nqfront,iel,nregi] = convex(neleg,npoig,ieleg,coorg,npfront,nqfront,ia,nregi);

% FILL IN ISIDE 
disp(' *** FILLING ISIDE FOR THE BACKGROUND GRID')
[iside,iloca,lwher,icone] = side(neleg,npoig,nsidg,ieleg,iside,lwher,lhowm,icone);

% FILL IN INTMEG
disp(' *** FILLING ELEMENT CONECTIVITY MATRIX')
nsidg = single(iloca);
intmeg = elecon(nsidg,ieleg,iside,intmeg);

xmax = single(-1.e+6);
xmin = single(1.e+6);
ymax = single(-1.e+6);
ymin = single(1.e+6);

% READ THE NEW MESH DATA
fid_17 = fopen(unit_17, 'r');
% get TEXT and NREG,NFN,NBCS 
temp = fscanf(fid_17,' %s %s %s',[3 1]);
temp = temp';
text = single(temp(:, :));

temp = fscanf(fid_17,' %f %f %f',[3 1]);
temp = temp';
nreg = single(temp(:, 1));
nfn  = single(temp(:, 2));
nbcs = single(temp(:, 3));


% SPECIFY COORDINATES OF FIXED NODES
 % get ieleg
 fid_17 = fopen(unit_17, 'r');
 fgets(fid_17);
 fgets(fid_17);
 fgets(fid_17);
 temp = fscanf(fid_17,'%f %f %f',[3 nfn]);
 temp = temp';
 
 ip    = single(temp(:,1));
 coorn = single(temp(:,2:3));
 coorn = coorn';
 
 % FIND OUT MAX AND MIN COORDINATES FOR PLOTTING
 for ip = 1:nfn
    xl=coorn(1,ip);
    yl=coorn(2,ip);

    if(xl < xmin)
        xmin=xl;
    end
    if(yl < ymin)
        ymin=yl;
    end
    if(xl > xmax)
        xmax=xl;
    end
    if(yl > ymax)
        ymax=yl;
    end
 end

inew = single(0);
coor0 = single([]);
% EXPAND BACKGROUND MESH TO INCLUDE ALL DOMAIN
coorg = expand(neleg,npoig,coorg,coor0,ieleg,intmeg,n1body,n2body,unit_7,fid_7,cv);
coorg = single(coorg);
iback = single(0);
ipoii= single(0);
ipoip = single(0);
node = single(0);

%  SPECIFY BOUNDARY SEGMENTS
disp(' *** GENERATING BOUNDARY NODES ')
nbou = single(0);
fgets(fid_17);
coor = single([]);
nbno = single([]);
nnn = single([]);
lbou = single([]);

lboud = single([]);
rcond = single([]);
igraph = single([]);

% ตัว๝ปลเฉพาะฝิจ
% t  = single(zeros(2,500));
% lx = single(zeros(500,1));
% c  = single(zeros(2,500));
% ti = single(zeros(500,1));
% ts = single(zeros(500,1));
% td = single(zeros(500,1));
% tx = single(zeros(500,1));
% xx = single(zeros(2,500));

t  = single([]);
lx = single([]);
c  = single([]);
ti = single([]);
ts = single([]);
td = single([]);
tx = single([]);
xx = single([]);

ni = single([]);
np = single([]);
tg = single([]);
tl = single([]);

%%%%%%%%%%% from findtin
xt = single([]);
arZ = single([]);

i1   =single([]);
i2   =single([]);
i3   =single([]);
ienr =single([]);
ind  =single([]);

tl = single([]);
srein = single([]);


amo  = single([]);
ttan = single([]);
dis  = single([]);
alp  = single([]);
anx  = single([]);
any  = single([]);

x = single([]);

ibcs = single([]);
nonf = single([]);
nonr = single([]);
iel = single([]);
ti = single([]);
ts = single([]);
td = single([]);
tx = single([]);
ti = single([]);
for ibs = 1:nbcs

    % GENERATE BOUNDARY NODES 

    [node, lboud, nbno, nbou, nnn, lbou, coor, dis, alp, anx, any, t, lx,c,ti,ts,td,tx,xx, ni, np, tg,arZ, i1,i2,i3,ienr,ind,tl, srein, xt,amo,ttan,x,rcond, ax, ay,xq,ilast] = interp(npoig,neleg,node,coorg,ieleg,intmeg,delta,coor,coorn,nbno,nnn,lcoor,lbou,nbou,ibs,ilast,ipoii,ipoip,lboud,rcond,igraph,unit_17,t, lx,c,ti,ts,td,tx,xx, ni, np, tg,xt,arZ,i1,i2,i3,ienr,ind,tl, srein, amo, dis,alp,anx,any,ttan,x);

end

nboun = node;

% TAKE CARE THE COMMON POINTS

% loop 30
for ibs = 1:nbcs

    % loop 34
    for i =1:2

        j = 1;

        if (i == 2)

            j = nnn(ibs);

        end

        kpoi = nbno(ibs,j);
        
        % loop 32
        for ibou = 1:nbou

            knew = lbou(2,ibou);
            ktry = lbou(1,ibou);
            
            if (ktry == kpoi)

                % goto 33
                break

            end
        end

        % come to 33
        nbno(ibs,j) = knew;
        
    end
end

% SPECIFY REGION BOUNDARY SEGMENTS
fid_17 = fopen(unit_17, 'r');
for i = 1:(3 + nbcs + 1 + 2*nbcs +1)
    fgetl(fid_17);
end
% loop 40
for ireg = 1:nreg
    
    % get mbcs
    temp = fscanf(fid_17,'%f %f',[1 2]);
    ip = single(temp(1));
    mbcs = single(temp(2));

    % get ibcs and mbcs
    temp = fscanf(fid_17,'%f %f %f %f %f',[mbcs 1]);
    temp = temp';
    ibcs = single(temp);
end




% ENQUIRE IF GENERATED NODES ARE REQUIRED

% SET UP STREC FOR THE ALREADY EXISTING NODES
disp(' *** INTERPOLATING FROM THE BACKGROUND GRID ')

% loop 1075
for i = 1:node

    x = coor(1,i);
    y = coor(2,i);
    
    inorm = 1;

    [ienr, arZ, i1, i2, i3,ilast] = findel(npoig,neleg,coorg,ieleg,intmeg,x,y,ilast,arZ,i1,i2,i3,ienr);

    [dis, alp, anx, any] = getvalue(npoig,i1,i2,i3,arZ,delta,dis,alp,anx,any,inorm);
    
    strec(1,i) = dis;
    strec(2,i) = alp;
    strec(3,i) = anx;
    strec(4,i) = any;
end

% SET UP INITIAL FRONTS

[nonf,nonr,nregi,npfront,nqfront] = setup(nreg,mbcs,ibcs,nbno,nnn,coor,node,npfront,nqfront,nonf,nregi,nonr);

% ประฝาศตัว๝ปลเฉพาะฝิจ

nelem =0;
toler = 0;

% TRIANGULATE REGIONS 
%sum(nelem,'all')
%sum(strec,'all')
igraph = 0 ;
% ENSURE SINGLE
nreg = single(nreg);
nonf = single(nonf);
npfront =single(npfront);
nqfront = single(nqfront);
nregi = single(nregi);
nonr = single(nonr);
iel = single(iel);
coor = single(coor);
nelem = single(nelem);
node = single(node);
neleg = single(neleg);
npoig = single(npoig);
ieleg = single(ieleg);
intmeg = single(intmeg);
coorg = single(coorg);
delta = single(delta);
strec = single(strec);
toler = single(toler);
ilast = single(ilast);
igraph = single(igraph);
[node, nelem, toler, iel,coor,fid_7] = triangle(nreg,nonf,npfront,nqfront,nregi,nonr,iel,coor,nelem,node,neleg,npoig,ieleg,intmeg,coorg,delta,strec,toler,ilast,igraph,unit_7,fid_7,cv,dis);

% place holder
% node = 1656
% nelem = 3150

% FIND OUT BOUNDARY POINTS
lcoor = boundar(node,nelem,lcoor,iel);


% FIND OUT REAL AND IDEAL NUMBER OF CONECTIVITIES
lcore = single(zeros(node,node));

lcore = conere(node,nelem,iel,lcore);

lcoid = single(zeros(node));

lcoid = coneid(node,nelem,iel,lcoor,lcoid,coor,strec);



% FILL IN ISIDE 
disp(' *** GENERATING BOUNDARY NODES ')

nside = nsidg;
[iside, iloca,lwher,icone] = side(nelem,node,nside,iel,iside,lwher,lhowm,icone);

nside = iloca;

fid_7 = fopen(unit_7, 'r');
for i =1:3
    fgets(fid_7);
end
%delta = delta';

counter = 0;
while true
    disp(' *** WHAT SHALL WE DO NOW ?')
    disp('     0 - QUIT')
    disp('     1 - PLOT THE MESH')
    disp('     2 - SMOOTH THE MESH')
    disp('     3 - SWAP DIAGONALS')
    disp('     4 - EAT 3 S ')
    disp('     5 - AREA CHECK/OUTPUT NO OF NODES AND ELEMENTS')
    disp('     6 - GET THE RE-START FILE')
    
    

    %temp = fscanf(fid_7,'%f',[1 1]);
    %temp = temp';
    %iwhat = temp(1);
    %iwhat = single(iwhat);
    %counter = counter + 1;
    intmat = iel(:,1:nelem);
    figure(1);
    z = coor(1,:).*0;
    trimesh(intmat',coor(1,:)',coor(2,:)',z,'EdgeColor','black');
    view(2),axis([xmin xmax ymin ymax]),axis equal, grid off, hold on
    %holder = input('Continue:? (Y/N) ', 's');
    iwhat = input(' Please select from the menu : \n');
    disp(iwhat)
    switch iwhat
        case 0
            disp('stop')
            break
        case 1
            intmat = iel(:,1:nelem);
            figure(1);
            z = coor(1,:).*0;
            trimesh(intmat',coor(1,:)',coor(2,:)',z,'EdgeColor','black');
            view(2),axis([xmin xmax ymin ymax]),axis equal, grid off, hold on
            %holder = input('Continue:? (Y/N) ', 's');
            continue
        case 2  
%             temp = fscanf(fid_7,'%f',[1 1]);
%             temp = temp';
%             nsmoo = temp(1);
%             nsmoo = single(nsmoo);
            nsmoo = input('ENTER NUMBER OF SMOOTHING LOOPS \n');
            coor = smooth(2,3,nelem,node,nsmoo,lcoor,lcore,iel,coor,coor0);
            continue
        case 3

            [iside, lcore, iel,nside] = swapdi(node,nelem,nside,iside,iel,lcore,lcoid,lcoor,coor,lwher,lhowm,icone);
            
            continue
        case 4
            lposi = [];

            [lcore,icone,coor,iside,lcoid,iel,nelem,nside] = eat3(node,nelem,nside,iel,iside,lcoor,lposi,lwher,lhowm,lcore,lcoid,icone,coor);
            
            continue
        case 5
            [coor] = areach(node,nelem,iel,coor);

            
            continue
        case 6
            output(npoig,neleg,coorg,ieleg,intmeg,unkng,node,nelem,coor,iel,rcond,nboun,lboud,iside,nside, toler,ilast,filnam,arZ,cv,delkmma,delkmmi,delsca)
            
            continue
        case 7
            continue
        otherwise
            disp('other value')
            continue
    end
end
return