function [erg] = ionex_read(pname)
% Read IONEX-files

url=fopen(strcat(pname));
while (1)
 txt=fgetl(url);
 chk=strfind(txt,'EPOCH OF FIRST MAP');
 if ~isempty(chk)
     break
 end
end
zwi.start_year=str2double(txt(1:6));
zwi.start_month=str2double(txt(7:12));
zwi.start_day=str2double(txt(13:18));
zwi.start_hour=str2double(txt(19:24));
zwi.start_minute=str2double(txt(25:30));
zwi.start_second=str2double(txt(31:36));

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'EPOCH OF LAST MAP');
 if ~isempty(chk)
     break
 end
end
zwi.end_year=str2double(txt(1:6));
zwi.end_month=str2double(txt(7:12));
zwi.end_day=str2double(txt(13:18));
zwi.end_hour=str2double(txt(19:24));
zwi.end_minute=str2double(txt(25:30));
zwi.end_second=str2double(txt(31:36));

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'INTERVAL');
 if ~isempty(chk)
     break
 end
end
zwi.interval=str2double(txt(1:6));


while (1)
 txt=fgetl(url);
 chk=strfind(txt,'# OF MAPS IN FILE');
 if ~isempty(chk)
     break
 end
end
zwi.number_of_maps=str2double(txt(1:6));

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'MAPPING FUNCTION');
 if ~isempty(chk)
     break
 end
end
zwi.mapping_function=txt(3:6);

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'ELEVATION CUTOFF');
 if ~isempty(chk)
     break
 end
end
zwi.cutoff=str2double(txt(1:8));

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'BASE RADIUS');
 if ~isempty(chk)
     break
 end
end
zwi.radius=str2double(txt(1:8));

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'MAP DIMENSION');
 if ~isempty(chk)
     break
 end
end
zwi.map_dimension=str2double(txt(1:6));

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'HGT1 / HGT2 / DHGT');
 if ~isempty(chk)
     break
 end
end

zwi.hgt1=str2double(txt(3:8));
zwi.hgt2=str2double(txt(10:15));
zwi.dhgt=str2double(txt(17:22));

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'LAT1 / LAT2 / DLAT');
 if ~isempty(chk)
     break
 end
end

zwi.lat1=str2double(txt(3:8));
zwi.lat2=str2double(txt(10:15));
zwi.dlat=str2double(txt(17:22));

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'LON1 / LON2 / DLON');
 if ~isempty(chk)
     break
 end
end

zwi.lon1=str2double(txt(3:8));
zwi.lon2=str2double(txt(10:15));
zwi.dlon=str2double(txt(17:22));

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'EXPONENT');
 if ~isempty(chk)
     break
 end
end

zwi.exponent=str2double(txt(1:6));

while (1)
 txt=fgetl(url);
 chk=strfind(txt,'END OF HEADER');
 if ~isempty(chk)
     break
 end
end

% starts reading VTEC maps
karten=[];
for i = 1:zwi.number_of_maps
    % find line START OF TEC MAP and read it "out"
    while (1)
      txt=fgetl(url);
      chk=strfind(txt,'START OF TEC MAP');
      if ~isempty(chk)
          break
      end
    end  
    % find line with the epoch and read it "out" 
    txt=fgetl(url);
    
    % read the epoch specific map
    karte=[];
    for j=1:length(zwi.lat1:zwi.dlat:zwi.lat2)
         while (1)
           txt=fgetl(url);
           chk=strfind(txt,'LAT/LON1/LON2/DLON/H');
           if ~isempty(chk)
             break 
           end
         end          
        lese_breite=[];
      % j ...   Index for latitude 
        L=length(zwi.lon1:zwi.dlon:zwi.lon2);        
        anz_zeilen=fix(L/16);
        if anz_zeilen>0
           for k=1:anz_zeilen
               txt=fgetl(url);
               lese_breite=[lese_breite str2double(txt(1:5))];
               lese_breite=[lese_breite str2double(txt(6:10))];
               lese_breite=[lese_breite str2double(txt(11:15))];
               lese_breite=[lese_breite str2double(txt(16:20))];
               lese_breite=[lese_breite str2double(txt(21:25))];
               lese_breite=[lese_breite str2double(txt(26:30))];
               lese_breite=[lese_breite str2double(txt(31:35))];
               lese_breite=[lese_breite str2double(txt(36:40))];
               lese_breite=[lese_breite str2double(txt(41:45))];
               lese_breite=[lese_breite str2double(txt(46:50))];
               lese_breite=[lese_breite str2double(txt(51:55))];
               lese_breite=[lese_breite str2double(txt(56:60))];
               lese_breite=[lese_breite str2double(txt(61:65))];
               lese_breite=[lese_breite str2double(txt(66:70))];
               lese_breite=[lese_breite str2double(txt(71:75))];
               lese_breite=[lese_breite str2double(txt(76:80))];
           end                      
        end % if anz_zeilen>0
        if mod(L,16)>0
            txt=fgetl(url);  
            for k=1:mod(L,16)              
              lese_breite=[lese_breite str2double(txt(1+(k-1)*5:5+(k-1)*5))];  
            end
        end % mod
         
     karte=[karte; lese_breite];    
    end  % Index j
  karten(:,:,i)=karte;
end
zwi.map=karten;
%fclose(url);

%url=fopen(strcat(pname,fname));


% starts reading RMS maps
karten=[];
for i = 1:zwi.number_of_maps
    %  find line START OF RMS MAP and read it "out"
    while (1)
      txt=fgetl(url);
      chk=strfind(txt,'START OF RMS MAP');
      chk2=strfind(txt,'END OF FILE');
      if ~isempty(chk)
          break
      end
      if ~isempty(chk2)
          break     
          txt
      end      
    end  
   if ~isempty(chk2)
       break;
   end
       
    % find line with the epoch and read it "out" 
    txt=fgetl(url);
    
    % read the epoch specific map
    karte=[];
    for j=1:length(zwi.lat1:zwi.dlat:zwi.lat2)
         while (1)
           txt=fgetl(url);
           chk=strfind(txt,'LAT/LON1/LON2/DLON/H');
           if ~isempty(chk)
             break 
           end
         end          
        lese_breite=[];
      % j ...   Index for Latitude
        L=length(zwi.lon1:zwi.dlon:zwi.lon2);        
        anz_zeilen=fix(L/16);
        if anz_zeilen>0
           for k=1:anz_zeilen
               txt=fgetl(url);
               lese_breite=[lese_breite str2double(txt(1:5))];
               lese_breite=[lese_breite str2double(txt(6:10))];
               lese_breite=[lese_breite str2double(txt(11:15))];
               lese_breite=[lese_breite str2double(txt(16:20))];
               lese_breite=[lese_breite str2double(txt(21:25))];
               lese_breite=[lese_breite str2double(txt(26:30))];
               lese_breite=[lese_breite str2double(txt(31:35))];
               lese_breite=[lese_breite str2double(txt(36:40))];
               lese_breite=[lese_breite str2double(txt(41:45))];
               lese_breite=[lese_breite str2double(txt(46:50))];
               lese_breite=[lese_breite str2double(txt(51:55))];
               lese_breite=[lese_breite str2double(txt(56:60))];
               lese_breite=[lese_breite str2double(txt(61:65))];
               lese_breite=[lese_breite str2double(txt(66:70))];
               lese_breite=[lese_breite str2double(txt(71:75))];
               lese_breite=[lese_breite str2double(txt(76:80))];
           end                      
        end % if anz_zeilen>0
        
        if mod(L,16)>0
            txt=fgetl(url);  
            for k=1:mod(L,16)              
              lese_breite=[lese_breite str2double(txt(1+(k-1)*5:5+(k-1)*5))];  
            end
        end % mod     
     karte=[karte; lese_breite];    
 end  % Index j        



  karten(:,:,i)=karte;
end % end von for i= 

% if isempty(chk2)    
  zwi.rms_map=karten;
% else
%   zwi.rms_map=99999;   
% end;
fclose(url);

erg=zwi;