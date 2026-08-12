% ************************************************************************
%   Description:
%	Manipulates the position of the satellite by changing one of the 
%   orbital elements (extending semi-major axis, increasing inclination) 
%   for orbit data from a FSO file.
%
%   Input:										
%     sources               sources struct
%     parameter             parameter struct
%     T2C_s                 transformation matrix to convert from TRS to CRS
%     numKepEle             number of keplerian element
%							 1: semi-major axis
%                            2: eccentricity
%                            3: inclination
%                            4: argument of perigee
%                            5: right ascension of ascending node
%                            6: argument of latitude
%
% 
%   Output:
%     sourcesChanged       changed sources struct
%     dKepEle			   change of orbital element
% 
%   External calls: 	
%
%       
%   Coded for VieVS: 
%   17 Dec 2024 by Helene Wolf
%
%   Revision: 
%
%
% ************************************************************************

function [sourceChanged, dKepEle] = manipulateFSOData(GM, source, T2C_s, numKepEle)
    
    sourceChanged = source;
    for k=1:length(source.x_crf)
        time = source.mjd(k);
        tosc = source.tosc;
        r = [source.x_crf(k), source.y_crf(k), source.z_crf(k)];
        v = [source.vx_crf(k), source.vy_crf(k), source.vz_crf(k)];

        [a,e,i,Omega,omega,t0] = xyzele(GM,time,r,v);

        dele=[0.01,0.000777,mas2rad(100),mas2rad(100),mas2rad(100),mas2rad(100)];
        u0 = getu0(GM, a,e,omega,t0,tosc);

        if numKepEle==1
            a = a+dele(1);
            t00 = gett0(GM,a,e,omega,tosc,u0);
            dKepEle = dele(1);
        elseif numKepEle==2
            e=e+dele(2);
            t00 = gett0(GM,a,e,omega,tosc,u0);
            dKepEle = dele(2);
        elseif numKepEle==3
            i=i+dele(3);
            t00 = gett0(GM,a,e,omega,tosc,u0);
            dKepEle = dele(3);
        elseif numKepEle==4
            Omega=Omega+dele(4);
            t00 = gett0(GM,a,e,omega,tosc,u0);
            dKepEle = dele(4);
        elseif numKepEle==5
            omega=omega+dele(5);
            t00 = gett0(GM,a,e,omega,tosc,u0);
            dKepEle = dele(5);
        elseif numKepEle==6
            t00 = gett0(GM,a,e,omega,tosc,u0+dele(6));
            dKepEle = dele(6);
        end
        
        [rn,vel1] = ephem(GM,a,e,i,Omega,omega,t00,time);
        sourceChanged.x_crf(k) = rn(1);
        sourceChanged.y_crf(k) = rn(2);
        sourceChanged.z_crf(k) = rn(3);
        sourceChanged.vx_crf(k) = vel1(1);
        sourceChanged.vy_crf(k) = vel1(2);
        sourceChanged.vz_crf(k) = vel1(3);
    end
    sourcesC.s(1)=sourceChanged;
    sourceChanged = getTRF_PosVelSatellites(sourcesC, T2C_s);
    sourceChanged = sourceChanged.s(1);
end