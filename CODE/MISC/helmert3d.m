function [tp,rc,ac,tr]=helmert3d(datum1,datum2,Type,WithOutScale,Approx,NameToSave)
% HERLMERT3D    overdetermined cartesian 3D similarity transformation ("Helmert-Transformation")
%
% [param, rotcent, accur, resid] = helmert3D(datum1,datum2,Type,DontUseScale,ApproxRot,SaveIt)
%
% Inputs:  datum1  n x 3 - matrix with coordinates in the origin datum (x y z)
%                  datum1 may also be a file name with ASCII data to be processed. No point IDs, only
%                  coordinates as if it was a matrix.
%
%          datum2  n x 3 - matrix with coordinates in the destination datum (x y z)
%                  datum2 may also be a file name with ASCII data to be processed. No point IDs, only
%                  coordinates as if it was a matrix.
%                  If some coordinates of datum2 are missing (e.g. a 2D point is used together with 
%                  3D points), NaN has to be used at the gap.
%                  If either datum1 and/or datum2 are ASCII files, make sure that they are of same
%                  length and contain corresponding points. there is no auto-assignment of points!
%
%            Type  is either '7p' for 7-parameter-transformation (Bursa-Wolf)
%                  or '10p' for 7-parameter-transformation with rotation center at centroid of
%                  datum1 (Molodensky-Badekas).
%                  Default is '7p' which is also choosen if Type is not a string (e.g. []).
%
%    DontUseScale  if this is not 0, do not calculate scale factor but set it to the inputted value
%                  Default: 0 (Use scale)
%
%       ApproxRot  1 x 3 approximate initial values for rotations. If rotation values are too big,
%                  the adjustment may fail if no or bad approximate values are given. Especially this
%                  is the case if rotation around Y-axis is close to multiples of pi/2. You might
%                  purport approximate values to overcome this weakness.
%                  Input as [ex ey ez] in radians.
%                  If omitted or left empty, default is [0 0 0].
%
%          SaveIt  string with the name to save the resulting parameters in Transformations.mat.
%                  Make sure only to use characters which are allowed in Matlab variable names
%                  (e.g. no spaces, umlaute, leading numbers etc.)
%                  If left out or empty, no storage is done.
%                  If the iteration is not converging, no storage is done and a warning is thrown.
%                  If the name to store already is existing, it is not overwritten automatically.
%                  To overwrite an existing name, add '#o' to the name, e.g. 'wgs84_to_local#o'.
%
% Outputs:  param  7 x 1 Parameter set of the 3D similarity transformation
%                      3 translations (x y z) in [Unit of datums]
%                      3 rotations (ex ey ez) in [rad]
%                      1 scale factor
%
%         rotcent  3 x 1 vector of rotation center in datum1 [x y z]
%
%           accur  7 x 1 accuracy of the parameters (or 6 x 1 if scale factor is set to be 1)
%
%           resid  n x 3 - matrix with the residuals datum2 - f(datum1,param)
%
% Used to calculate datum transformation parameters e.g. for global change of reference ellipsoid
% when at least 3 identical points in both datum systems are known.
% 
% Systems need to be right-handed, i.e. [x y z]. Makes use of Tait-Bryan 
% angle order XYZ for passive intrinsic rotations (see explanation document
% for this toolbox for formulas)- it is EPSG:9607.
%
% ATTENTION: Please be aware of the approximate values-problem mentioned above. In the scope of
%            datum transformation, this is no problem as rotation angles are sufficiently small.
%            This function is trying to use an affine approach to determine proper approximations
%            if more than 3 ID points are given and throws a warning if results become unsecure.
%            (Results may be ambiguous, while either solution is of equal accuracy, though. There
%            may be different sets of rotations leading to nearly the same result.) If no approximate
%            values are given with only 3 ID points, a warning may be thrown if big rotations are
%            probably. But don't rely on that; please check the results for plausibility in any case.
%            If no approximate information is known, but more than 3 ID points are given, you might
%            use an affine transformation (helmertaffine3d) instead. This is tried by the function
%            automatically. 
%
% Parameters can be used with d3trafo.m
% This function needs helmertaffine3d.m
% 04/14/09 Peter Wasmeier - Technische Universität München
% p.wasmeier@bv.tum.de
% 04/17/15 Andrea Pinna - University degli Studi di Cagliari
% andrea.pinna@unica.it
% Corrected a small typo in the derivatives thanks to Marten from Neubrandenburg
 
%% Argument checking and defaults
if nargin<6
    NameToSave=[];
else
    NameToSave=strtrim(NameToSave);
    if strcmp(NameToSave,'#o')
        NameToSave=[];
    end
end
if nargin<5 || isempty(Approx)
    Approx=[0 0 0];
elseif numel(Approx)~=3
    error('ApproxRot needs to be a 3-element vector.')
else
    Approx=Approx(:)';
end
if nargin<4 || isempty(WithOutScale)
    WithOutScale=0;
end
if nargin<3 || ~ischar(Type) || isempty(Type)
    Type='7p'; 
end
% Load input file if specified
if ischar(datum1)
    datum1=load(datum1);
end
if ischar(datum2)
    datum2=load(datum2);
end
if (size(datum1,1)==3)&&(size(datum1,2)~=3)
    datum1=datum1'; 
end
if (size(datum2,1)==3)&&(size(datum2,2)~=3)
    datum2=datum2'; 
end
s1=size(datum1);
s2=size(datum2);
if any(s1~=s2)
    error('The datum sets are not of equal size')
elseif any([s1(2) s2(2)]~=[3 3])
    error('At least one of the datum sets is not 3D')
elseif any([s1(1) s2(1)]<3)
    error('At least 3 points in each datum are necessary for calculating')
elseif any(isnan(datum1))
    error('NaN are not allowed in the origin datum; complete 3D coordinates are necessary there.')
elseif prod(s2)-sum(sum(isnan(datum2))) < 7
    error('At least 7 coordinates in the destination datum are necessary')
end
switch Type
    case '7p'
        rc=[0 0 0];
    case '10p'
        rc=mean(datum1);
    otherwise
        error ('Transformation type needs to be ''7p'' or ''10p''.')
end
%% Approximation values
naeh=[[mean(datum2(~isnan(datum2(:,1)),1)) mean(datum2(~isnan(datum2(:,2)),2)) mean(datum2(~isnan(datum2(:,3)),3))]-mean(datum1) Approx 1];
if all(Approx==[0 0 0]) && s1(1)>3
    try
        x0=helmertaffine3d(datum1,datum2);
        s=(sqrt(x0(4)^2+x0(5)^2+x0(6)^2)+sqrt(x0(7)^2+x0(8)^2+x0(9)^2)+sqrt(x0(10)^2+x0(11)^2+x0(12)^2))/3;
        if abs(x0(11))<1e-6 && abs(x0(12))<1e-6
            if x0(10)<0
                ey=3*pi/2;
            else
                ey=pi/2;
            end
            warning('Helmert3D:Ambiguous_rotations','Y-rotation is close to a multiple of pi/2. X- and Z-rotation therefore cannot be approximated.')
            ex=0;
            ez=0; 
        else
            ex=atan2(-x0(11),x0(12));
            ey=atan2(x0(10),sqrt((x0(4))^2+(x0(7))^2));
            ez=atan2(-x0(7),x0(4));
        end
    catch
        ex=0;
        ey=0;
        ez=0;
        s=1;
    end 
    naeh=[0 0 0 ex ey ez s];
end
if WithOutScale
    naeh(7)=WithOutScale;
end
WertA=[1e-5 1e-8];
zaehl=0;
x0=naeh(1);
y0=naeh(2);
z0=naeh(3);
ex=naeh(4);
ey=naeh(5);
ez=naeh(6);
m=naeh(7);
tp=[x0 y0 z0 ex ey ez m];
% This will only be used, if point accuracies were given - not used yet
% Qbb=speye(3*s1(1));
%% Adjustment
SuppressSingularMatrixWarning=0;
while(1)
    
    % Compute sine and cosine of angles - speeds up computation
	cosex = cos(ex);
	cosey = cos(ey);
	cosez = cos(ez);
	sinex = sin(ex);
	siney = sin(ey);
	sinez = sin(ez); 
                
    A=zeros(3*s1(1),7);
    w=zeros(3*s1(1),1);
    A(1:3:end,1) = -1;
    A(2:3:end,2) = -1;
    A(3:3:end,3) = -1;
    
    A(1:3:end,4)=-m*((cosex*siney*cosez-sinex*sinez)*(datum1(:,2)-rc(2))+(sinex*siney*cosez+cosex*sinez)*(datum1(:,3)-rc(3)));
    A(1:3:end,5)=-m*((-siney*cosez)*(datum1(:,1)-rc(1))+(sinex*cosey*cosez)*(datum1(:,2)-rc(2))+(-cosex*cosey*cosez)*(datum1(:,3)-rc(3)));
    A(1:3:end,6)=-m*((-cosey*sinez)*(datum1(:,1)-rc(1))+(-sinex*siney*sinez+cosex*cosez)*(datum1(:,2)-rc(2))+(+cosex*siney*sinez+sinex*cosez)*(datum1(:,3)-rc(3)));
    A(1:3:end,7)=-((cosey*cosez)*(datum1(:,1)-rc(1))+(sinex*siney*cosez+cosex*sinez)*(datum1(:,2)-rc(2))+(-cosex*siney*cosez+sinex*sinez)*(datum1(:,3)-rc(3)));
    
    A(2:3:end,4)=-m*((-cosex*siney*sinez-sinex*cosez)*(datum1(:,2)-rc(2))+(-sinex*siney*sinez+cosex*cosez)*(datum1(:,3)-rc(3)));
    A(2:3:end,5)=-m*((siney*sinez)*(datum1(:,1)-rc(1))+(-sinex*cosey*sinez)*(datum1(:,2)-rc(2))+(cosex*cosey*sinez)*(datum1(:,3)-rc(3)));
    A(2:3:end,6)=-m*((-cosey*cosez)*(datum1(:,1)-rc(1))+(-sinex*siney*cosez-cosex*sinez)*(datum1(:,2)-rc(2))+(cosex*siney*cosez-sinex*sinez)*(datum1(:,3)-rc(3)));
    A(2:3:end,7)=-((-cosey*sinez)*(datum1(:,1)-rc(1))+(-sinex*siney*sinez+cosex*cosez)*(datum1(:,2)-rc(2))+(cosex*siney*sinez+sinex*cosez)*(datum1(:,3)-rc(3)));
    
    A(3:3:end,4)=-m*((-cosex*cosey)*(datum1(:,2)-rc(2))+(-sinex*cosey)*(datum1(:,3)-rc(3)));
    A(3:3:end,5)=-m*((cosey)*(datum1(:,1)-rc(1))+(-sinex*(-siney))*(datum1(:,2)-rc(2))+(-cosex*siney)*(datum1(:,3)-rc(3)));
    A(3:3:end,6)=0;
    A(3:3:end,7)=-((siney)*(datum1(:,1)-rc(1))+(-sinex*cosey)*(datum1(:,2)-rc(2))+(cosex*cosey)*(datum1(:,3)-rc(3)));
   
    w(1:3:end)=-rc(1)+datum2(:,1)-x0-m*((cosey*cosez)*(datum1(:,1)-rc(1))+(sinex*siney*cosez+cosex*sinez)*(datum1(:,2)-rc(2))+(-cosex*siney*cosez+sinex*sinez)*(datum1(:,3)-rc(3)));
    w(2:3:end)=-rc(2)+datum2(:,2)-y0-m*((-cosey*sinez)*(datum1(:,1)-rc(1))+(-sinex*siney*sinez+cosex*cosez)*(datum1(:,2)-rc(2))+(cosex*siney*sinez+sinex*cosez)*(datum1(:,3)-rc(3)));
    w(3:3:end)=-rc(3)+datum2(:,3)-z0-m*((siney)*(datum1(:,1)-rc(1))+(-sinex*cosey)*(datum1(:,2)-rc(2))+(cosex*cosey)*(datum1(:,3)-rc(3)));
        
    if WithOutScale
        A=A(:,1:6);
    end
    
    nans=find(isnan(datum2'));
    A(nans,:)=[];
    w(nans)=[]; 
    
    lastwarn('');
    warning('off','MATLAB:nearlySingularMatrix');
    w=-1*w;
    r=size(A,1)-size(A,2);
    %Pbb=inv(Qbb);
    Pbb = speye(size(A,1));   % Accuracy of points is not used yet
    APbbA = A'*Pbb*A;
    deltax=APbbA\(A'*Pbb*w);
    v=A*deltax-w;
    sig0p=sqrt((v'*Pbb*v)/r);
    Kxxda=(sig0p^2)*inv(APbbA);
    ac=sqrt(diag(Kxxda));
    warning('on','MATLAB:nearlySingularMatrix');
    
    % Test for warning
    [WarnMessages,WarnID]=lastwarn;
    if SuppressSingularMatrixWarning == 0 && strcmp(WarnID,'MATLAB:nearlySingularMatrix')
        WarnMessages=[WarnMessages '\nThis does not necessarily mean that results are inaccurate, but sometimes they may be.'];
        if strcmp(Type,'7p')
            WarnMessages=[WarnMessages '\nConsider changing the transformation type to 10 parameters to overcome this.'];
        else
            WarnMessages=[WarnMessages '\nPlease check whether your input data is sufficient for computing (also geometrical point distribution).'];
        end
        warning('Helmert3D:NearlySingularMatrix',WarnMessages);
        SuppressSingularMatrixWarning = 1;
    end
    testv=sqrt((deltax(1)^2+deltax(2)^2+deltax(3)^2)/3);
    testd=sqrt((deltax(4)^2+deltax(5)^2+deltax(6)^2)/3);
    zaehl=zaehl+1;
    x0=x0+deltax(1);
    y0=y0+deltax(2);
    z0=z0+deltax(3);
    ex=ex+deltax(4);
    ey=ey+deltax(5);
    ez=ez+deltax(6);
    if ~WithOutScale && (m+deltax(7))>1e-15     % This condition is to prevent numerical problems with m-->0
        m=m+deltax(7);
    end
    tp=[x0 y0 z0 ex ey ez m]';
    if abs(testv) < WertA(1) && abs(testd) < WertA(2)
        break;
    elseif zaehl>1000
        warning('Helmert3D:Too_many_iterations','Calculation not converging after 1000 iterations. I am aborting. Results may be inaccurate.')
        break;
    end
end
if any (abs(tp(4:6))>2*pi)
    warning('Helmert3D:Unsufficient_approximation_values','Rotation angles seem to be big. A better approximation is regarded. Results will be inaccurate.')
end
%% Transformation residuals
idz=zeros(s1);
idz(:,2)=rc(2)+tp(2)+tp(7)*((-cos(tp(5))*sin(tp(6)))*(datum1(:,1)-rc(1))+(-sin(tp(4))*sin(tp(5))*sin(tp(6))+cos(tp(4))*cos(tp(6)))*(datum1(:,2)-rc(2))+(cos(tp(4))*sin(tp(5))*sin(tp(6))+sin(tp(4))*cos(tp(6)))*(datum1(:,3)-rc(3)));
idz(:,1)=rc(1)+tp(1)+tp(7)*((cos(tp(5))*cos(tp(6)))*(datum1(:,1)-rc(1))+(sin(tp(4))*sin(tp(5))*cos(tp(6))+cos(tp(4))*sin(tp(6)))*(datum1(:,2)-rc(2))+(-cos(tp(4))*sin(tp(5))*cos(tp(6))+sin(tp(4))*sin(tp(6)))*(datum1(:,3)-rc(3)));
idz(:,3)=rc(3)+tp(3)+tp(7)*((sin(tp(5)))*(datum1(:,1)-rc(1))+(-sin(tp(4))*cos(tp(5)))*(datum1(:,2)-rc(2))+(cos(tp(4))*cos(tp(5)))*(datum1(:,3)-rc(3)));
tr=datum2-idz;
%% Store data
if ~isempty(NameToSave)
    load Transformations;
    if zaehl>1000
        warning('Helmert3D:Results_too_inaccurate_to_save','Results may be inaccurate and do not get stored.')
    elseif exist(NameToSave,'var') && length(NameToSave)>=2 && ~strcmp(NameToSave(end-1:end),'#o')
        warning('Helmert3D:Parameter_already_exists',['Parameter set ',NameToSave,' already exists and therefore is not stored.'])
    else
        if strcmp(NameToSave(end-1:end),'#o')
            NameToSave=NameToSave(1:end-2);
        end
        if any(regexp(NameToSave,'\W')) || any(regexp(NameToSave(1),'\d'))
            warning('Helmert3D:Parameter_name_invalid',['Name ',NameToSave,' contains invalid characters and therefore is not stored.'])
        else
            eval([NameToSave,'=[',num2str(tp'),' ',num2str(rc),'];']);
            save('Transformations.mat',NameToSave,'-append');
        end
    end
end