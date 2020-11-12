

% Design matrix A(20) for the baseline dependent clock offsets

    
  function [A_bdco,H_bdco,Ph_bdco,ebsl] = abdo_clk(scan,antenna,opt,n_observ,parameter)
      
    % all baselines from the scan   
    bslst1=[]; bslst2=[];
    for i=1:length(scan)
        bslst1 = [bslst1 scan(i).obs.i1];
        bslst2 = [bslst2 scan(i).obs.i2];
    end


    na = length(antenna); r=[];
    if opt.bdco_fromOPT == 0 % don't read the list from the OPT file, but create it now
        % list of all posible baseline combinations
        bl1=[]; bl2=[];
        for i=1:na-1
            bl0=[];
            bl0(1:na-i,1) = i;
            bl1 = [bl1; bl0];

            bl0=[];
            bl0(1:na-i,1) = [i+1:na];
            bl2 = [bl2; bl0];
        end
        ebsl = [bl1 bl2]; 
        
        
        % compute number of observations per baseline in the "ebsl" list
        nobs_bsl = zeros(size(ebsl,1),1);
        for i=1:size(ebsl,1)
           nobs_bsl(i) = length (intersect ( find(bslst1 == ebsl(i,1)) , find(bslst2 == ebsl(i,2)) )) +...
                         length (intersect ( find(bslst1 == ebsl(i,2)) , find(bslst2 == ebsl(i,1)) ));
        end

        % take the reference station from the OPT file "REFERENCE CLOCK"
        rstat = find([opt.stat.ref] == 1);
        idref = [find(bl1 == rstat ); find(bl2 == rstat )];         
        
        % check if there are baselines with less than opt.bdco_minobs observations to the
        % reference station (in GUI opt.bdco_minobs = 5)
        idfewrefobs = find(nobs_bsl(idref) < opt.bdco_minobs);
        
        % if there are baselines with less than opt.bdco_minobs
        % observations, search for another reference station. It is the
        % station with the highest number of observations
        badrefstat=[]; statwithobs=1;
        while ~isempty(idfewrefobs) && statwithobs == 1
           badrefstat = [badrefstat; rstat; ebsl(idref(idfewrefobs),1); ebsl(idref(idfewrefobs),2)];
           badrefstat = unique(badrefstat);
           
           % possible reference stations
            posiblestat = setdiff([1:na],badrefstat);
            if isempty(posiblestat)
                statwithobs=0;
%                error('The algorithm for BCs did not find any station with observations at all baselines.')
           end
           
           noall(1:na,1)=0;           
           noall(posiblestat,1) = [antenna(posiblestat).numobs];
           
           rstat =[];
           [~, rstat] = (max(noall)); % new reference station
           
           idref=[]; idfewrefobs = [];
           idref = [find(bl1 == rstat ); find(bl2 == rstat )];
           idfewrefobs = find(nobs_bsl(idref) < opt.bdco_minobs); % check the new reference station
        end
        
        koef = 0;
        if statwithobs == 0
            % as reference take the station with max observations
            rstat =[]; idref=[]; idfewrefobs = [];
            [~, rstat] = max([antenna.numobs]);
            
            % fix station with the minimum of observations
            minrstat =[];
            [~, minrstat] = min([antenna.numobs]);
            
            % find the baselines with the reference stat and the stat with
            % minimum observations - will be fixed
            idrefI = [find(bl1 == rstat ); find(bl2 == rstat )];    
            idminstat = [find(bl1 == minrstat ); find(bl2 == minrstat )];
            
            % find the station with the minumum of baselines - will be
            % fixed
            for i = 1:na
               nrbas(i)= length(nonzeros(nobs_bsl([find(bl1 == i ); find(bl2 == i )])));                
            end
            [~, minnrbas] = min(nrbas);
            idminnrbas = [find(bl1 == minnrbas ); find(bl2 == minnrbas )];
            
            idminnrbas2=[];
            if minrstat == minnrbas
               % find the station with the second lowest nr of baselines 
                nrbas(minnrbas)=100;
                [~, minnrbas2] = min(nrbas);
                idminnrbas2 = [find(bl1 == minnrbas2 ); find(bl2 == minnrbas2 )];
            end
            
            % increase the minimal number of observations per baseline for 10
            koef = 10;

            idref = [idrefI; idminstat; idminnrbas; idminnrbas2];
            idref = unique(idref);
        end
        

        % delete baselines from the "ebsl" list with the chosen stations
        ebsl(idref,:) = [];
        nobs_bsl(idref,:) = [];

        % delete baselines from the "ebsl" list with less than parameter.lsmopt.bdco_minobs + koef observations
        ebsl(find(nobs_bsl < (opt.bdco_minobs + koef)),:) = []; % 
    
    else
        for i=1:na
            id=[];
            id=strcmp(antenna(i).name, {parameter.opt.options.bdco_est.sta1});
            ebsla(id,1)=i;
            id=[];
            id=strcmp(antenna(i).name, {parameter.opt.options.bdco_est.sta2});
            ebsla(id,2)=i;
        end
        % delete baselines if the station was not found in antenna.mat 
        r=[];  [r,~]=find(ebsla==0);
        ebsla(r,:)=[]; 
        
        % switch the order of stations in the baseline if necessary
        idchange = find(ebsla(:,1)>ebsla(:,2)); 
        ebsl = ebsla;
        ebsl(idchange,1) = ebsla(idchange,2);
        ebsl(idchange,2) = ebsla(idchange,1);
    end

    
    A_bdco=[];H_bdco=[];Ph_bdco=[];
    for ib = 1:size(ebsl,1) % loop over selected baselines (ebsl) where the bas-dep clock offset will be estimated
        ibsl = []; ibsl1 = []; ibsl2 = [];

        
        ibsl1 = (intersect ( find(bslst1 == ebsl(ib,1)) , find(bslst2 == ebsl(ib,2)) ));
        ibsl2 = (intersect ( find(bslst1 == ebsl(ib,2)) , find(bslst2 == ebsl(ib,1)) ));
       
        ibsl = [ibsl1 ibsl2];
        
        Abslclk = [];
        Abslclk = zeros(n_observ,1);
        Abslclk(ibsl,1) = 1;

        % Concatenating 
        A_bdco = horzcat(A_bdco, Abslclk);

        H_bdco(ib) = 0;
        Ph_bdco(ib) = 0;

    end
    
    fprintf('Baseline-dependent clock offset: %1.0f\n', size(ebsl,1))
    for k=1:size(ebsl,1)
        fprintf('%8s   %8s \n', antenna(ebsl(k,1)).name, antenna(ebsl(k,2)).name)
    end
    if ~isempty(r)
    fprintf('Following baselines (given in OPT file) not found:\n')
        for k=1:length(r)
            fprintf('%8s   %8s \n', parameter.opt.options.bdco_est(r(k)).sta1, parameter.opt.options.bdco_est(r(k)).sta2)
        end
    end

