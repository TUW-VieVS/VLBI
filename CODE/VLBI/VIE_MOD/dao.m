function [d,grn,gre] = dao (stat8,az,el)
%
% # MGR-FORMAT V 1.0   GSFC DAO + GMAO Mean Site Gradients 
% *************************************************************
% *
% * Mean site gradients computed from GSFC Data Assimilation
% * Office (DAO) model for met data from 1990-95. 
% * Values for missing and new stations were taken from mean
% * site gradients computed from the GSFC GMAO (Global Modeling
% * and Assimilation Office) GEOS-5 numerical weather model data
% * for data from 2011.0-2013.0.
% * Atmospheric * gradient delay is modeled as 
% *      tau = m_grad(el,az) * [GN*cos(az)+GE*sin(az)], 
% * where el and az are the elevation and azimuth of the the 
% * observation and the gradient mapping function is m_grad.
% * The gradient vector has east and north components
% * GE and GN. Here, the delay units are mm.
% *
% * References:
% *
% * MacMillan, D.S., Atmospheric gradients from very
% *    long baseline interferometry observations, Geophys. Res.
% *    Lett., 22, 1041-1044, 1995.
% *
% * MacMillan, D.S. and C. Ma, Atmospheric gradients and the 
% *    VLBI terrestrial and celestial reference frames, Geophys.
% *    Res. Lett., 24, 453-456, 1997.
% *
% *
% *  pet  2001.09.25  Added KOGANEI, MOURA, TATEYAMA: copied values from KASHIMA
% *  pet  2002.05.20  Added AIRA, CHICHI10, GIFU3:    copied values from KASHIMA
% *  pet  2002.05.29  Added TIGOCONC: copied values from SANTIA12
% *  pet  2003.03.06  Added CTVASBAY, CTVASTJ   copied from ALGOPARK
% *  pet  2003.06.27  Added VERAMZSW 
% *  pet  2003.07.21  Added GIFU11, GBT-VLBA
% *  pet  2005.03.25  Added METSAHOV ( copy of METSHOVI )
% *  pet  2005.03.25  Added SVETLOE  ( copy of METSHOVI )
% *  pet  2005.12.27  Added TIDBIN64 ( copy of DSS15 )
% *  pet  2006.01.19  Added ZELENCHK ( copy of CRIMEA )
% *  pet  2006.02.09  Added VERAIRIK 
% *  pet  2006.02.09  Added VERAISGK 
% *  pet  2006.02.09  Added VERAOGSW
% *  pet  2007.03.23  Added DSS13
% *  pet  2007.03.23  Added BADARY
% *  DG   2009.04.01  Added GOLDMARS and MADRID64, corrected TIDBIN64 (copy of DSS45)
% *  DG   2010.04.15  Added YEBES40M and SINTOTU3.
% *  DM   2010.10.21  Added HOBART12 (copied HOBART26)
% *  DM   2014.02.04  Added UCHINOUR (copied MIZUGSI)
% *  DM   2014.02.04  Added KUNMING (copied SESHAN25)
% *  DM   2014.02.04  Added HART15M (copied HARTRAO)
% *  DM   2014.02.04  Added KATH12M,YARRA12M,WARK12M (from GMAO_2011_13)
% *  DM   2014.02.04  Added KOGANEI,KOGANEI3,MIURA,TATEYAMA,
% *                   TIDBIN64,KAINAN,KANOZAN,NOBEYAMA,SUWON,
% *                   TOMAKO11 (from GMAO_2011_13)
% *  DM   2014.02.04  Updated with values from GMAO_2011_13:
% *                   AIRA,CHICH10,CTVASBAY,CTVASTJ,MIZNAO10,
% *                   MIZUSGSI,SVETLOE,TIGOCONC
% *  DM   2015.02.13  Added TIANMA65 (copied SESHAN25)
% *  DM   2015.02.18  Added WETTZ13N (copied WETTZELL)
% *  DM   2015.02.18  Added WETTZ13S (copied WETTZELL)
% *  DM   2015.02.18  Added VLA-N8 (copied VLA)
% *  DM   2015.02.19  Added SEJONG (copied MIYAZAKI)
% *  DM   2015.03.04  Added RAEGYEB (copied YEBES)
% *
% *************************************************************
%
% input data
% ----------
% stat8: 8-character name of the VLBI station
% az: azimuth from north in radians
% el: elevation angle in radians
% 
% output data
% -----------
% d: asymmetric delay in meter
% grn: north gradient in mm
% gre: east gradient in mm
% 
% Johannes Boehm, 2011 Feb 8
% update: Hana Krasna, 2014 Feb 10; GSFC update file in Feb 2014
%
% 16 May 2014 Hana; Typo by SYOWA corrected
% 04 Jan 2016 David; Updated file (new file on server)
% 
daotxt = [...
 'AIRA        -.56E+00    -.25E-01 ';
 'ALGOPARK    -.46E+00    -.15E-01 ';
 'AUSTINTX    -.25E+00     .46E-01 ';
 'AZORES      -.11E+00    -.71E-01 ';
 'BADARY      -.26E+00    -.58E-01 ';
 'BERMUDA     -.26E+00     .42E-02 ';
 'BLKBUTTE    -.14E+00     .82E-01 ';
 'BLOOMIND    -.40E+00    -.20E-03 ';
 'BR-VLBA     -.17E+00     .31E-01 ';
 'BREST       -.27E+00    -.55E-02 ';
 'CARNUSTY    -.35E+00     .17E-01 ';
 'CARROLGA    -.28E+00     .21E-01 ';
 'CHICHI10    -.25E+00    -.24E+00 ';
 'CHLBOLTN    -.28E+00    -.34E-02 ';
 'CRIMEA      -.16E+00     .58E-03 ';
 'CTVASBAY    -.49E+00     .73E-01 ';
 'CTVASTJ     -.46E+00     .12E+00 ';
 'DAITO       -.18E+00    -.21E-01 ';
 'DEADMANL    -.88E-01     .77E-01 ';
 'DSS13       -.86E-01     .67E-01 ';
 'DSS15       -.86E-01     .67E-01 ';
 'DSS45        .24E+00     .73E-02 ';
 'DSS65       -.14E+00     .41E-02 ';
 'DSS65A      -.14E+00     .41E-02 ';
 'EFLSBERG    -.30E+00     .10E-01 ';
 'ELY         -.16E+00     .15E-01 ';
 'FD-VLBA     -.24E+00    -.19E-01 ';
 'FLAGSTAF    -.15E+00     .70E-01 ';
 'FORT ORD    -.11E+00     .42E-01 ';
 'FORTLEZA     .18E-01    -.47E-01 ';
 'FORTORDS    -.11E+00     .40E-01 ';
 'FTD 7900    -.24E+00    -.19E-01 ';
 'GBT-VLBA    -.42E+00    -.21E-01 ';
 'GGAO7108    -.40E+00    -.36E-01 ';
 'GIFU11      -.50E+00     .32E-01 ';
 'GIFU3       -.50E+00     .32E-01 ';
 'GILCREEK    -.18E+00     .34E-01 ';
 'GOLDECHO    -.23E+00     .14E+00 ';
 'GOLDMARS    -.86E-01     .67E-01 ';
 'GOLDVENU    -.86E-01     .68E-01 ';
 'GORF7102    -.40E+00    -.36E-01 ';
 'GRASSE      -.80E-01    -.32E-01 ';
 'HALEAKAL    -.13E+00    -.79E-01 ';
 'HART15M      .10E+00    -.11E-02 ';
 'HARTRAO      .10E+00    -.11E-02 ';
 'HATCREEK    -.13E+00     .66E-01 ';
 'HAYSTACK    -.43E+00    -.24E-01 ';
 'HN-VLBA     -.43E+00    -.22E-01 ';
 'HOBART12     .29E+00     .81E-02 ';
 'HOBART26     .29E+00     .81E-02 ';
 'HOFN        -.25E+00     .34E-01 ';
 'HOHENFRG    -.31E+00    -.26E-01 ';
 'HOHNBERG    -.19E+00    -.28E-01 ';
 'HRAS 085    -.24E+00    -.19E-01 ';
 'JPL MV1     -.95E-01     .54E-01 ';
 'KAINAN      -.57E+00     .18E+00 ';
 'KANOZAN     -.64E+00    -.14E+00 ';
 'KARLBURG    -.30E+00    -.28E-01 ';
 'KASHIM11    -.50E+00     .32E-01 ';
 'KASHIM34    -.50E+00     .32E-01 ';
 'KASHIMA     -.50E+00     .32E-01 ';
 'KATH12M      .29E+00    -.22E-01 ';
 'KAUAI       -.19E+00    -.22E-01 ';
 'KIRSBERG    -.29E+00    -.36E-01 ';
 'KODIAK      -.16E+00     .47E-01 ';
 'KOGANEI     -.43E+00    -.18E+00 ';
 'KOGANEI3    -.43E+00    -.18E+00 ';
 'KOKEE       -.19E+00    -.22E-01 ';
 'KP-VLBA     -.19E+00     .13E+00 ';
 'KUNMING     -.61E+00    -.11E+00 ';
 'KWAJAL26    -.29E+00    -.91E-01 ';
 'LA-VLBA     -.17E+00     .13E-01 ';
 'LEONRDOK    -.31E+00     .28E-01 ';
 'MADRID64    -.14E+00     .41E-02 ';
 'MAMMOTHL    -.84E-01     .66E-01 ';
 'MARCUS      -.15E+00     .12E-02 ';
 'MARPOINT    -.40E+00    -.44E-01 ';
 'MATERA      -.11E+00    -.34E-01 ';
 'MCD 7850    -.24E+00    -.18E-01 ';
 'MEDICINA    -.93E-01    -.32E-01 ';
 'METSAHOV    -.29E+00    -.50E-02 ';
 'METSHOVI    -.29E+00    -.50E-02 ';
 'MIAMI20     -.18E+00     .25E-01 ';
 'MILESMON    -.32E+00    -.81E-01 ';
 'MIURA       -.65E+00    -.12E+00 ';
 'MIYAZAKI    -.53E+00     .21E-01 ';
 'MIZNAO10    -.49E+00    -.92E-01 ';
 'MIZUSGSI    -.49E+00    -.92E-01 ';
 'MK-VLBA     -.84E-01    -.88E-01 ';
 'MOJ 7288    -.86E-01     .67E-01 ';
 'MOJAVE12    -.86E-01     .67E-01 ';
 'MON PEAK    -.14E+00     .89E-01 ';
 'MV2ONSLA    -.29E+00    -.13E-01 ';
 'NL-VLBA     -.37E+00    -.35E-01 ';
 'NOBEY 6M    -.51E+00     .63E-01 ';
 'NOBEYAMA    -.51E+00     .63E-01 ';
 'NOME        -.11E+00     .37E-01 ';
 'NOTO        -.15E+00    -.83E-01 ';
 'NRAO 140    -.42E+00    -.21E-01 ';
 'NRAO20      -.42E+00    -.21E-01 ';
 'NRAO85 1    -.42E+00    -.21E-01 ';
 'NRAO85 3    -.42E+00    -.22E-01 ';
 'NYALES20    -.19E+00     .22E-02 ';
 'OCOTILLO    -.14E+00     .90E-01 ';
 'OHIGGINS     .42E+00    -.30E-01 ';
 'ONSALA60    -.29E+00    -.13E-01 ';
 'ONSALA85    -.29E+00    -.13E-01 ';
 'OV-VLBA     -.90E-01     .70E-01 ';
 'OVR 7853    -.90E-01     .70E-01 ';
 'OVRO 130    -.90E-01     .70E-01 ';
 'PARKES       .25E+00     .14E-01 ';
 'PBLOSSOM    -.91E-01     .59E-01 ';
 'PENTICTN    -.17E+00     .27E-01 ';
 'PIETOWN     -.18E+00     .49E-01 ';
 'PINFLATS    -.14E+00     .83E-01 ';
 'PLATTVIL    -.24E+00    -.44E-01 ';
 'PRESIDIO    -.12E+00     .64E-01 ';
 'PT REYES    -.99E-01    -.15E-01 ';
 'PVERDES     -.15E+00     .48E-01 ';
 'QUINCY      -.95E-01     .73E-01 ';
 'RAEGYEB     -.15E+00     .40E-02 ';
 'RICHMOND    -.18E+00     .25E-01 ';
 'ROBLED32    -.14E+00     .41E-02 ';
 'SAGARA      -.51E+00     .48E-01 ';
 'SANPAULA    -.11E+00     .57E-01 ';
 'SANTIA12     .24E+00     .19E+00 ';
 'SC-VLBA     -.89E-01    -.12E+00 ';
 'SEATTLE1    -.24E+00     .35E-01 ';
 'SEJONG      -.53E+00     .21E-01 ';
 'SESHAN25    -.61E+00    -.11E+00 ';
 'SEST         .14E+00     .14E+00 ';
 'SHANGHAI    -.62E+00    -.11E+00 ';
 'SINTOTU     -.43E+00     .23E-01 ';
 'SINTOTU3    -.43E+00     .23E-01 ';
 'SNDPOINT    -.22E+00     .61E-01 ';
 'SOURDOGH    -.14E+00     .50E-01 ';
 'SUWON       -.59E+00     .79E-01 ';
 'SVETLOE     -.29E+00     .55E-02 ';
 'SYOWA        .00E+00     .00E+00 '; 
 'TATEYAMA    -.64E+00    -.14E+00 ';
 'TIANMA65    -.61E+00    -.11E+00 ';
 'TIDBIN64     .26E+00     .83E-01 ';
 'TIGOCONC     .14E+00     .23E+00 ';
 'TIGOWTZL    -.26E+00    -.44E-01 ';
 'TITIJIMA    -.23E+00    -.14E-02 ';
 'TOMAKO11    -.41E+00    -.10E+00 ';
 'TOULOUSE    -.17E+00    -.40E-01 ';
 'TROMSONO    -.30E+00     .51E-02 ';
 'TRYSILNO    -.24E+00    -.73E-02 ';
 'TSUKU3      -.55E+00     .32E-01 ';
 'TSUKUB32    -.55E+00     .32E-01 ';
 'TSUKUBA     -.55E+00     .32E-01 ';
 'UCHINOUR    -.49E+00    -.92E-01 ';
 'URUMQI      -.26E+00    -.58E-01 ';
 'USSURISK    -.37E+00    -.53E-01 ';
 'USUDA64     -.55E+00     .65E-01 ';
 'VERAIRIK    -.57E+00     .25E-01 ';
 'VERAISGK    -.57E+00     .25E-01 ';
 'VERAMZSW    -.57E+00     .25E-01 ';
 'VERAOGSW    -.57E+00     .25E-01 ';
 'VERNAL      -.20E+00     .22E-01 ';
 'VICTORIA    -.18E+00     .36E-01 ';
 'VLA         -.18E+00     .50E-01 ';
 'VLA-N8      -.18E+00     .50E-01 ';
 'VNDNBERG    -.14E+00     .40E-02 ';
 'WARK12M      .42E+00     .71E-01 ';
 'WESTFORD    -.43E+00    -.24E-01 ';
 'WETTZ13N    -.26E+00    -.44E-01 ';
 'WETTZ13S    -.26E+00    -.44E-01 ';
 'WETTZELL    -.26E+00    -.44E-01 ';
 'WHTHORSE    -.15E+00     .40E-01 ';
 'YAKATAGA    -.86E-01     .72E-01 ';
 'YARRA12M     .45E+00     .13E+00 ';
 'YEBES       -.15E+00     .40E-02 ';
 'YEBES40M    -.15E+00     .40E-02 ';
 'YELLOWKN    -.29E+00    -.68E-01 ';
 'YLOW7296    -.29E+00    -.68E-01 ';
 'YUMA        -.18E+00     .78E-01 ';
 'ZELENCHK    -.16E+00     .58E-03 '];
                                     
                                     
                                     
[n,m] = size(daotxt);                
                                     
grnlist = zeros(n,1);               
grelist = zeros(n,1);
statlist = '';
 
for i = 1:n
    statlist(i,:) = daotxt(i,1:8);
    grnlist(i) = str2double(daotxt(i,11:20));
    grelist(i) = str2double(daotxt(i,23:32));
end
 
% north and east gradients
grn = 0.d0;
gre = 0.d0;
for i = 1:n
    if strcmp(stat8,statlist(i,1:8))
        grn = grnlist(i);
        gre = grelist(i);
    end
end
 
% calculation of the asymmetric delay in m (Chen and Herring 1997)
d = 1/(sin(el)*tan(el)+0.0032)*(grn*cos(az)+gre*sin(az))/1000; 
 
