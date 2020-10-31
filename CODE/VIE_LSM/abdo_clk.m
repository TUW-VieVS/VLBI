

% Design matrix A(20) for the baseline dependent clock offsets

    
  function [A_bdco,H_bdco,Ph_bdco,ebsl] = abdo_clk(scan,antenna,opt,n_observ,parameter)
      
    % all baselines from the scan
    u = struct2cell(scan);
    u6 = u(6,:,:);
    bslst1=[]; bslst2=[];
    for i=1:length(scan)
        bslst1 = [bslst1 [u6{i}.i1]];
        bslst2 = [bslst2 [u6{i}.i2]];
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

        % delete baselines from the "ebsl" list with the reference station
        rstat = find([opt.stat.ref] == 1);
        idref = [find(bl1 == rstat ); find(bl2 == rstat )]; 
        ebsl(idref,:) = [];


        % compute number of observations per baseline in the "ebsl" list
        nobs_bsl = zeros(size(ebsl,1),1);
        for i=1:size(ebsl,1)
           nobs_bsl(i) = length (intersect ( find(bslst1 == ebsl(i,1)) , find(bslst2 == ebsl(i,2)) ));
        end
        % delete baselines from the "ebsl" list with less than parameter.lsmopt.bdco_minobs observations
        ebsl(find(nobs_bsl < opt.bdco_minobs),:) = []; % 
    
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
        f1 = []; f2 = []; ibsl = [];

        f1 = find(bslst1==ebsl(ib,1));
        f2 = find(bslst2==ebsl(ib,2));

        ibsl = intersect(f1, f2);

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
        fprintf('%8s - %8s \n', antenna(ebsl(k,1)).name, antenna(ebsl(k,2)).name)
    end
    if ~isempty(r)
    fprintf('Following baselines (given in OPT file) not found:\n')
        for k=1:length(r)
            fprintf('%8s - %8s \n', parameter.opt.options.bdco_est(r(k)).sta1, parameter.opt.options.bdco_est(r(k)).sta2)
        end
    end

