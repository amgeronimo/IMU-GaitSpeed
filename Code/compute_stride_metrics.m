function [] = compute_stride_metrics(metadata,redcap_path,proc_path,RCdata,RCind,record_id_vec,Vbl)

subj = Vbl.subj;
accHeader = Vbl.accHeader;
gyrHeader = Vbl.gyrHeader;
pathtodata = Vbl.pathtodata;

for si = 1:length(subj)
    
    kk=1;
    clear oSL oDU oWS oNM oDA oSA oSS oSE oSLthresh_valid
    
    kp = strcmpi({metadata.name},subj{si});
    metadatatmp = metadata(kp);
    
    
    
    %Order dates
    serialdates = zeros(length(metadatatmp),1);
    for i = 1:length(metadatatmp)
        serialdates(i) = datenum(str2num(metadatatmp(i).date(1:4)),str2num(metadatatmp(i).date(5:6)),str2num(metadatatmp(i).date(7:8)));
    end
    % find min
    [~,loc] = min(serialdates);
    if loc ~= 1
        new_start = loc;
        metatmp = metadatatmp(loc:end);
        metatmp2 = metadatatmp(1:loc-1);
        metadatatmp = [metatmp,metatmp2]; % new metadata
    end
    
    clear ssn snm snc
    for k = 1:length(metadatatmp)
        ssn{k} = metadatatmp(k).session;
        snm{k} = metadatatmp(k).name;      % name for particular session
        snc{k} = metadatatmp(k).name(regexpi(metadatatmp(k).name,'[^SC]'));      % just the study code
    end
    
    [unm,ind] = unique(ssn,'rows','stable'); % Omit duplicate results
    [ind0,ind2,ind3] = unique(snc); % was using snm but now scn to account for cases when one subject's study code is changed
    sctr = 0;
    oNM_old = '';
    for i = ind'
        stepFol = [redcap_path.bigOrg,filesep,redcap_path.indiv,filesep,...
            metadatatmp(i).SCname,filesep,metadatatmp(i).SCname,proc_path.step,filesep];
        %     clear stridestats
        try
            SSdata = load([stepFol,'M2parameters_imu_',metadatatmp(i).session,'.mat']);
            stridestats = SSdata.stridestats;
        catch
            disp(['M2 results for ' metadatatmp(i).session ' missing'])
            stridestats = [];
        end
        
        oNM{kk} = metadatatmp(i).SCname;
        oSS{kk} = [metadatatmp(i).SCname,'-',metadatatmp(i).date(5:6),'/',metadatatmp(i).date(7:8)];
        oSA{kk} = metadatatmp(i).assist;
        
        if strcmp(oNM(kk),oNM_old)
            sctr = sctr+.05;
        else
            sctr = 0;
        end
        
        oDA{kk} = metadatatmp(i).date;
        
        %If more space is needed betweeen the boxplots, of different
        %individuals, ind3(i) can be multiplied by a factor greater than 1.
        oSE{kk} = ind3(i)+sctr;
        oNM_old = oNM{kk};
        
        if isempty(stridestats)
            oSL{kk} = [];
            oDU{kk} = [];
            oWS{kk} = [];
            oSLthresh_valid{kk} = [];
            
        else
            %         tmp = zeros(1,length(stridestats.SL));
            %         tmp(stridestats.SL<.1)=1;
            %         eh = error_ellipse_COMPILED(stridestats.SL(~tmp)',stridestats.DU(~tmp)',0);
            %         tmp2 = eh.outlier;
            %         tmp(tmp==0) = 2*tmp2;
            %         tmpcmp = [0 0 1;1 0 0;1 0 1];
            %         tmpcmp = tmpcmp(tmp+1,:);
            
            %         oSL{kk} = stridestats.SL(~tmp)';
            %         oDU{kk} = stridestats.DU(~tmp)';
            oSL{kk} = stridestats.SL';
            oDU{kk} = stridestats.DU';
            oWS{kk} = oSL{kk}./oDU{kk};
            %         oSLthresh_valid{kk} = tmp;
        end
        kk=kk+1;
    end
    
        
    %% Save and Export Data
    % this saves stride metricsc for an individual subject
    GSA_data = [];
    GSA_data.filesep = filesep;
    GSA_data.Headers = [accHeader;gyrHeader];
    GSA_data.pathtodata = pathtodata;
    GSA_data.proc_path = proc_path;
    GSA_data.RCdata = RCdata;
    GSA_data.RCind = RCind;
    GSA_data.record_id_vec = record_id_vec;
    GSA_data.redcap_path = redcap_path;
    GSA_data.oSS = oSS; 
    GSA_data.oSA = oSA;
    GSA_data.oDA = oDA;
    GSA_data.oNM = oNM; 
    GSA_data.oSE = oSE; 
    GSA_data.oSL = oSL;
    GSA_data.oDU = oDU;
    GSA_data.oWS = oWS;
    
    saveFol = [redcap_path.bigOrg,filesep,redcap_path.indiv,filesep,metadatatmp(1).SCname filesep];
    save([saveFol 'GSA_data_',char(datetime(clock,'Format','yyyy-MM-dd')), '_', subj{si},'.mat'],'GSA_data');
end