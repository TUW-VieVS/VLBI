% Comments:
% 1. fopen: avoid 't'
% 2. Output is little-endian
% 3. Output does not use FORTRAN style record markers.

% Binary file containing all information about source coordinates 
% used by Chris Jacobs

% created for VieVS by Hana Krasna


function newcrf_binary(globsol,paths)
% paths=pathGS % separate function
% clc
% load('globsol.mat')


% add the FORTRAN style record markers
FORTRANsrm = 1;


% c     Binary file containing all information about source coordinates
% c
% c     Four records:
% c
% c       1. Header record
% c          I*4  number of sources
% c          I*4  flag for covariance to follow (1 = yes)
% c          C*60 file name of SRIF output file

numSou = size(globsol.source.refname.IERS,1);


flagCov = 1; % 0/1
inFile = blanks(60);
inFile1 = [paths.L2]; % your description of the solution


inFile(1:length(inFile1)) = inFile1;

% paths.path_out='../OUT/GLOB/';
fileID = fopen([paths.path_out 'CRF/' paths.out '/crf_catalogue_binary_' paths.L2 '.catb'],'w');

%%
% [filename, permission, machineformat, encoding] = fopen(fileID)

if FORTRANsrm ==1; fwrite(fileID,68,'int32');end
fwrite(fileID,numSou,'int32');
fwrite(fileID,flagCov,'int32');
fwrite(fileID,inFile,'char');
if FORTRANsrm ==1; fwrite(fileID,68,'int32');end

% c       2. Source record (in RA order)
RA_apr=globsol.source.apriori_rade(:,1);  %[rad]
dRA = globsol.source.drade(:,1); %[ms]
RA_est=RA_apr + dRA/1000/3600/12*pi; %[rad]

De_apr=globsol.source.apriori_rade(:,2); %[rad]
dDe = globsol.source.drade(:,2); %[mas]
De_est=De_apr + dDe/1000/3600/180*pi; %[rad]

mRADe = globsol.source.sigma_rade; %[ms, mas]
mRA = mRADe(:,1)/1000/3600/12*pi; %[rad]
mDe = mRADe(:,2)/1000/3600/180*pi; %[rad]

Q = globsol.Q;
dradecol = globsol.source.dradecol;
souactiv = globsol.source.souactiv;

if FORTRANsrm ==1; fwrite(fileID,92*numSou,'int32');end
for k=1:numSou
    
    % c          C*12 source name
	src = blanks(12);
    src1 = globsol.source.refname.IVS(k,:);
    src(1:length(src1)) = src1;
    
    fwrite(fileID,src,'char');
    
    % c          R*8  right ascension, declination (radians)
    fwrite(fileID,RA_est(k),'float64');
    fwrite(fileID,De_est(k),'float64');
    
% c          R*8  right ascension, declination sigmas (radians)
%     mRA(k)
    fwrite(fileID,mRA(k),'float64');
    fwrite(fileID,mDe(k),'float64');
    
    % c          R*8  RA-dec correlation coefficient
    correl = Q(dradecol(k,1),dradecol(k,2)) / sqrt(Q(dradecol(k,1),dradecol(k,1)) * Q(dradecol(k,2),dradecol(k,2)) );
    fwrite(fileID,correl,'float64');
    
    % c          R*8  mean observation epoch (years)
    % c          R*8  earliest observation time (years)
    % c          R*8  last observation time (years)
    
    idObsSou=find(souactiv(k,:)>0);
    soumjdall = souactiv(end,idObsSou);
    
    first_last_mean_O_MJD = [min(soumjdall) max(soumjdall) (min(soumjdall)+max(soumjdall))/2] ;
    
    [yyDoySecod]=mjd2yydoysecod(first_last_mean_O_MJD(3));
    mean_epoch=yyDoySecod(1)+yyDoySecod(2)/365.25;
    fwrite(fileID,mean_epoch,'float64');
    
    [yyDoySecod]=mjd2yydoysecod(first_last_mean_O_MJD(1));
    timlo=yyDoySecod(1)+yyDoySecod(2)/365.25;
    fwrite(fileID,timlo,'float64');

    [yyDoySecod]=mjd2yydoysecod(first_last_mean_O_MJD(2));
    timhi=yyDoySecod(1)+yyDoySecod(2)/365.25;
    fwrite(fileID,timhi,'float64');

    BL3.mean_epoch(k) = mean_epoch;
    
    % c          I*4  number of observing sessions
    Nexp = length(idObsSou);
    fwrite(fileID,Nexp,'int32');

    % c          I*4  number of delay
    Nobs = sum(souactiv(k,idObsSou));
    fwrite(fileID,Nobs,'int32');
    
    % c          I*4  rate, phase observations
    fwrite(fileID,0,'int32');
    fwrite(fileID,0,'int32');
end
if FORTRANsrm ==1; fwrite(fileID,92*numSou,'int32');end


% c       3. Parameter record
% c          I*4  total number of estimated source parameters
% c          C*24 parameter names
% c          R*8  value of parameter, sigma
% c          R*8  reference time for parameter (years)

ntot = numSou *2; %RA and De
% number of parameters
if FORTRANsrm ==1; fwrite(fileID,4+ntot*48,'int32');end
fwrite(fileID,ntot,'int32');

% Loop over the parameter names
for k=1:numSou
    src = blanks(12);
    src1 = globsol.source.refname.IVS(k,:);
    src(1:length(src1)) = src1;

    p1 = ['RIGHT ASCEN.' src];
    p2 = ['DECLINATION.' src];

    fwrite(fileID,p1,'char');
    fwrite(fileID,p2,'char');
end

% Loop over the triplet of parameters
for k=1:numSou
    fwrite(fileID,RA_est(k),'float64');
    fwrite(fileID,mRA(k),'float64');
    fwrite(fileID,BL3.mean_epoch(k),'float64');
    fwrite(fileID,De_est(k),'float64');
    fwrite(fileID,mDe(k),'float64');
    fwrite(fileID,BL3.mean_epoch(k),'float64');
end
if FORTRANsrm ==1; fwrite(fileID,4+ntot*48,'int32');end

% c       4. Covariance record
% c          R*8  covariance matrix, upper triangle


if flagCov ==1
    if FORTRANsrm ==1; fwrite(fileID,(ntot*(ntot+1)/2)*8,'int32');end
    Q_sinex = globsol.Q;
    col_ra = globsol.source.dradecol(:,1);
    col_de = globsol.source.dradecol(:,2);
    c_sou = [col_ra; col_de];

    Q_sinex(c_sou,c_sou) = Q_sinex(c_sou,c_sou)./(1000*3600*180/pi)^2; %mas^2 --> rad^2

    tmpMatSou=globsol.source.dradecol';
    % make one vector out of it...
    tmpMat=[tmpMatSou(:)];
    % ... and delete zeros
    tmpMat(tmpMat==0)=[];

    % reorder the normal equation matrix
    % preallocate
    Q=zeros(length(tmpMat), length(tmpMat));

    for col=1:size(Q,1)
        for sou=1:numSou
            Q((sou-1)*2+1, col)=Q_sinex(col_ra(sou), tmpMat(col));
            Q((sou-1)*2+2, col)=Q_sinex(col_de(sou), tmpMat(col));
        end
    end
    
    % write normal equation matrix
    for col=1:size(Q,1)
        for row=1:col
%             Q(col,row);
            fwrite(fileID,Q(col,row),'float64');
        end
    end
    
    if FORTRANsrm ==1; fwrite(fileID,(ntot*(ntot+1)/2)*8,'int32');end
    
end



fclose(fileID);


% 
% 
% %% Read the file
% 
% FORTRANsrm=1; % FORTRAN style record marks
% 
% fileID = fopen('test1.catb');
% 
% if FORTRANsrm ==1; fread(fileID,1,'int32');end
% nso=fread(fileID,1,'int32');
% flagCov=fread(fileID,1,'int32');
% fread(fileID,60,'uint8=>char')';
% if FORTRANsrm ==1; fread(fileID,1,'int32');end
% 
% 
% if FORTRANsrm ==1; fread(fileID,1,'int32');end
% for i=1:nso
%     souname(i,:)=fread(fileID,12,'uint8=>char')';
%     RA(i)=fread(fileID,1,'float64');
%     De(i)=fread(fileID,1,'float64');
%     mRA(i)=fread(fileID,1,'float64');
%     mDe(i)=fread(fileID,1,'float64');
%     CorCoef_B3(i)=fread(fileID,1,'float64');
%     mean_epoch(i)=fread(fileID,1,'float64');
%     timlo(i)=fread(fileID,1,'float64');
%     timhi(i)=fread(fileID,1,'float64');
%     fread(fileID,1,'int32');
%     fread(fileID,1,'int32');
%     fread(fileID,1,'int32');
%     fread(fileID,1,'int32');
%         
% end
% if FORTRANsrm ==1; fread(fileID,1,'int32');end
% 
% 
% 
% 
% 
% % c       3. Parameter record
% % c          I*4  total number of estimated source parameters
% % c          C*24 parameter names
% % c          R*8  value of parameter, sigma
% % c          R*8  reference time for parameter (years)
% 
% if FORTRANsrm ==1; fread(fileID,1,'int32');end
% % Loop over the parameter names
% fread(fileID,1,'int32');
% for k=1:nso
%     fread(fileID,24,'uint8=>char')';
%     fread(fileID,24,'uint8=>char')';
% end
% 
% % Loop over the triplet of parameters
% for k=1:nso
%     
%     fread(fileID,1,'float64');
%     fread(fileID,1,'float64');
%     fread(fileID,1,'float64');
%     fread(fileID,1,'float64');
%     fread(fileID,1,'float64');
%     fread(fileID,1,'float64');
% 
% end
% if FORTRANsrm ==1; fread(fileID,1,'int32');end
% 
% % c       4. Covariance record
% % c          R*8  covariance matrix, upper triangle
%               
% if flagCov==1
%     if FORTRANsrm ==1; fread(fileID,1,'int32');end
%     
%     ntot=nso*2;
%     
%     for i=1:ntot
%         for j=1:i
%             C(j,i)=fread(fileID,1,'float64');
%         end
%     end
%     
%     if FORTRANsrm ==1; fread(fileID,1,'int32');end
% end
% 
% fclose(fileID);
% 
% for i=1:nso
%     correl_fromB4(i) = C(2*i-1,2*i) / sqrt(C(2*i-1,2*i-1) * C(2*i,2*i) );
% end