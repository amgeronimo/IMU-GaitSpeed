function [text] = paper_results(redcap_path,proc_path,Vbl)

m2_path = [redcap_path.bigOrg '\' filesep '\' redcap_path.indiv];
pathtohomerec = '..\HomeRecordings';


%% Load results from the main Redcap files
Vbl.pathtodata = ['..' filesep 'HomeRecordings' filesep];
Vbl.pathtoredcap = {['..' filesep 'Redcap_data' filesep 'redcap_ControlData.csv'];
    ['..' filesep 'Redcap_data' filesep 'redcap_PatientData.csv']};
[RC,RCind,record_id_vec] = redcap_input(redcap_path,proc_path,Vbl);

RCp = RC(cellfun(@(x) ~isempty(x),regexp({RC.study_code},'SC')));
RCc = RC(cellfun(@(x) ~isempty(x),regexp({RC.study_code},'^C')));

%Isolate subjects who have at least one 10mw task performed
keepsub = find(~cellfun(@isempty,{RC.gam_walk}));

%Find valid 10mw tests for these subjects
clear GA
kj=1;
for i = keepsub
    keepses = find(cellfun(@(x) sum(isempty(x))==0,RC(i).gam_walk));
    for j = keepses
        GA(kj).code = [RC(i).study_code];
        GA(kj).date = RC(i).gam_ga_date{j};
        GA(kj).tmw = [RC(i).gam_walk{j}{:}];
        GA(kj).tmw_speed = 6./GA(kj).tmw;
        GA(kj).tug = [RC(i).gam_tug{j}{:}];
        kj=kj+1;
    end
end

text.demo = [num2str(length(RCp)) ' patients (' num2str(nansum(cell2mat([RCp.e_sex]))),...
    ' male) and ' num2str(length(RCc)) ' controls (',...
    num2str(nansum(cell2mat([RCc.e_sex]))) ' male) were enrolled in the study.',...
    ];

text.gam = [num2str(length(GA)) ' of those enrolled completed ' num2str(length({GA.date})),...
    ' Gait and Motor Assessments'];

% % extract FRS scores for all clinic dates for all patients
% clear frs_all frs_gm frs_gm_sum frs_walk
% for i = 1:length(RCp)
%     for j = 1:length(RCp(i).c_alsfrs)
%         frs_all{i,j} = RCp(i).c_alsfrs{j};
%         frs_walk(i,j) = RCp(i).c_alsfrs{j}(8);
%         frs_gm{i,j} = RCp(i).c_alsfrs{j}(7:9);
%         frs_gm_sum(i,j) = sum(cell2mat(RCp(i).c_alsfrs{j}(7:9)));
%     end
% end

% % extract strength scores for all patients, watch out for NaNs in RC
% clear str_hipR str_hipL str_kflexR str_kflexL str_kextR str_kextL
% clear str_footR str_footL str_sum str_score str_avg str_avg_score
% for i = 1:length(RCp)
%     for j = 1:length(RCp(i).c_str_date)
%         str_hipR(i,j) = RCp(i).c_str_hip_r{j};
%         str_hipL(i,j) = RCp(i).c_str_hip_l{j};
%         str_kflexR(i,j) = RCp(i).c_str_k_flex_r{j};
%         str_kflexL(i,j) = RCp(i).c_str_k_flex_l{j};
%         str_kextR(i,j) = RCp(i).c_str_k_ext_r{j};
%         str_kextL(i,j) = RCp(i).c_str_k_ext_l{j};
%         str_footR(i,j) = RCp(i).c_str_foot_r{j};
%         str_footL(i,j) = RCp(i).c_str_foot_l{j};
%         % sum
%         str_sum(i,j) = sum([str_hipR(i,j),str_hipL(i,j),str_kflexR(i,j),str_kflexL(i,j),str_kextR(i,j),str_kextL(i,j),str_footR(i,j),str_footL(i,j)],'omitnan');
%         str_score(i,j) = str_sum(i,j)/40;
%         str_avg(i,j) = mean([str_hipR(i,j),str_hipL(i,j),str_kflexR(i,j),str_kflexL(i,j),str_kextR(i,j),str_kextL(i,j),str_footR(i,j),str_footL(i,j)],'omitnan');
%         str_avg_score(i,j) = str_avg(i,j)/5;
%     end
% end
%
% % days from first - clinic visits
% % FRS and strength
% clear frs_datenum frs_daysfromfirst
% clear str_datenum str_daysfromfirst
% for i = 1:length(RCp)
%     frs_datenum{i} = cellfun(@(x) datenum(x,'mm/dd/yy'),RCp(i).c_alsfrs_date);
%     frs_daysfromfirst{i} = frs_datenum{i}-frs_datenum{i}(1);
%     good_str = find(~ismissing(string(RCp(i).c_str_date)));
%     str_datenum{i} = cellfun(@(x) datenum(x,'mm/dd/yy'),RCp(i).c_str_date(good_str));
%     str_daysfromfirst{i} = str_datenum{i}-str_datenum{i}(1);
% end


% Demographics by group
% sex: 0 = female, 1 = male
% patients
ageP = zeros(length(RCp),1);
sexP = zeros(length(RCp),1);
htP = zeros(length(RCp),1);
wtP = zeros(length(RCp),1);
for i = 1:length(RCp)
    ageP(i) = RCp(i).e_age;
    sexP(i) = cell2mat(RCp(i).e_sex);
    htP(i) = cell2mat(RCp(i).c_height(1));
    wtP(i) = cell2mat(RCp(i).c_weight(1));
end
% backfill height for SC7
htP(3) = cell2mat(RCp(3).c_height(2));
% stats
% mean_ageP = mean(ageP);
med_ageP = median(ageP);
min_ageP = min(ageP);
max_ageP = max(ageP);
% std_ageP = std(ageP);
femalesP = length(sexP(sexP==0));
malesP = length(sexP(sexP==1));
% mean_wtP = mean(wtP);
% std_wtP = std(wtP);
med_wtP = median(wtP);
% mean_htP = mean(htP);
% std_htP = std(htP);
med_htP = median(htP);
% table
% mean_statsP = [mean_ageP,min_ageP,max_ageP,std_ageP,mean_wtP,std_wtP,mean_htP,std_htP,femalesP,malesP];
% mean_namesP = {'Mean Age','Min Age','Max Age','SD Age','Mean Wt','SD Wt','Mean Ht','SD Ht','#F','#M'};
% mean_tblP = array2table(mean_statsP,'VariableNames',mean_namesP)
med_statsP = [med_ageP,min_ageP,max_ageP,med_wtP,med_htP,femalesP,malesP];
med_namesP = {'Median Age','Min Age','Max Age','Median Wt','Median Ht','#F','#M'};
med_tblP = array2table(med_statsP,'VariableNames',med_namesP)

% controls
ageC = zeros(length(RCc),1);
sexC = zeros(length(RCc),1);
for i = 1:length(RCc)
    ageC(i) = RCc(i).e_age;
    sexC(i) = cell2mat(RCc(i).e_sex);
end
% stats
mean_ageC = mean(ageC);
min_ageC = min(ageC);
max_ageC = max(ageC);
std_ageC = std(ageC);
med_ageC = median(ageC);
femalesC = length(sexC(sexC==0));
malesC = length(sexC(sexC==1));
% table
statsC = [mean_ageC,min_ageC,max_ageC,med_ageC,std_ageC,femalesC,malesC];
namesC = {'Mean Age','Min Age','Max Age','Median Age','SD Age','#F','#M'};
tbleC = array2table(statsC,'VariableNames',namesC)

%% Stats on the home recordings
filedir = dir(pathtohomerec);

%Find files containing an underscore
tmp2 = regexp({filedir(:).name},'_');
goodind = ~cellfun(@isempty,tmp2);
%Find files ending in csv
tmp2 = regexp({filedir(:).name},'.csv');
goodind2 = ~cellfun(@isempty,tmp2);

badfiles =  {filedir(~goodind|~goodind2).name};
datafiles = {filedir(goodind&goodind2).name};

D.name = regexp(datafiles,'^[a-zA-Z0-9]+','match','once');
tmpt = regexp(D.name,'[a-zA-Z]*','match');
tmpn = regexp(D.name,'\d*','match');
D.SCname = cellfun(@(x,y) sprintf('%s%02d',upper(x),str2num(y)),[tmpt{:}],[tmpn{:}],'UniformOutput',false);
D.dates = regexp(datafiles,'[0-9]+-[0-9]+','match','once');
D.datenum = cellfun(@(x) datenum(x,'yyyymmdd-HHMMSS'),D.dates);
D.u_name = unique(D.SCname);
D.u_date = unique(D.datenum);

clear daydiff dates recrangedays numrecs dayfromfirst recweek change_assist ind_change change_name assistdata
for i = 1:length(D.u_name)
    ind = strcmp(D.SCname,D.u_name{i});

    [datesort sortind] = sort(D.datenum(ind));
    daydiff{i} = diff(datesort);

    %duplicates
    dups = daydiff{i}==0; % was <0.1 but failed for subjects who actually sent 2 recordings on the same day (C04,C05)
    dates{i} = datesort([true ~dups]); % will save this version
    datetmp = datesort([true ~dups]); % use this for analysis below
    daydiff{i}(dups) = [];

    recrangedays(i) = range(datetmp);
    numrecs(i) = length(unique(datetmp));
    dayfromfirst{i} = datetmp-datetmp(1);
    dayfromfirst{i} = unique(dayfromfirst{i});
    dayfromfirst_round{i} = round(dayfromfirst{i});
    recweek{i} = round(dayfromfirst{i}/7);

    % assist device
    assist_files = datafiles(ind);
    assist_files = assist_files(sortind);
    assist_unique = assist_files([true ~dups]);
    for j = 1:length(assist_unique)
        tmpund = cell2mat(strfind(assist_unique(j),'_'));
        tmpchar = char(assist_unique(j));
        assistdata{i}(j) = str2double(tmpchar(tmpund(1)+1));
    end
    % find subjects who changed assist device
    tmpdiff = assistdata{i}-circshift(assistdata{i},-1);
    if ~isempty(find(tmpdiff))
        tmpname = D.u_name{i};
        change_assist(i) = str2double(tmpname(regexpi(char(D.u_name{i}),'[^SC0]')));
        ind_change(i) = i;
        change_name{i} = tmpname;
    end
end

% remove zeros from change_assist
change_assist = change_assist(change_assist~=0);
ind_change = ind_change(ind_change~=0);

% take specific assist for the ones who changed
clear assist_sub
for i = 1:length(change_assist)
    assist_sub{i} = [change_assist(i);sort(assistdata{ind_change(i)})']; % because subject 11 does not have files in order, assuming assit number increases
end
% change assist device based on iphone or android
% for sc7 and sc8, no changes
% for sc11, decrease by 1
assist_sub{3}(2:end) = assist_sub{3}(2:end)-1;
assist_sub{3}(assist_sub{3}<0) = 0;

alldaydiff = cat(2,daydiff{:});
missedweeks = cellfun(@(x) sum(~ismember(min(x):max(x),x)),recweek);

patient_recs = numrecs(7:end);
controls_recs = numrecs(1:6);  %USE A SUBSET OF THESE
patients_range = recrangedays(7:end);
controls_range = recrangedays(1:6);
patients_missed = missedweeks(7:end);
controls_missed = missedweeks(1:6);
patients_alldays = alldaydiff(7:end);
controls_alldays = alldaydiff(1:6);
text.recs = [num2str(length(numrecs)) ' subjects completed a median of ',...
    num2str(median(numrecs)) ' (' num2str(min(numrecs)) '-' num2str(max(numrecs)),...
    ') recordings over ' sprintf('%.0f',median(recrangedays)) ' (',...
    sprintf('%.0f',min(recrangedays)),...
    '-' sprintf('%.0f',max(recrangedays)) ') days. ',...
    'All but ' num2str(sum(missedweeks>0)) ' subject(s) submitted at least',...
    ' one recording per week while in the study; in total only ',...
    sprintf('%.1f',100*sum(alldaydiff>=7)/length(alldaydiff)) '% of recordings',...
    ' were sent 7 or more days after the previous recording'];
text.patients = [num2str(length(patient_recs)) ' patients completed a median of ',...
    num2str(median(patient_recs)) ' (' num2str(min(patient_recs)) '-' num2str(max(patient_recs)),...
    ') recordings over ' sprintf('%.0f',median(patients_range)) ' (',...
    sprintf('%.0f',min(patients_range)),...
    '-' sprintf('%.0f',max(patients_range)) ') days. ',...
    'All but ' num2str(sum(patients_missed>0)) ' subject(s) submitted at least',...
    ' one recording per week while in the study; in total only ',...
    sprintf('%.1f',100*sum(patients_alldays>=7)/length(patients_alldays)) '% of recordings',...
    ' were sent 7 or more days after the previous recording'];
text.controls = [num2str(length(controls_recs)) ' controls completed a median of ',...
    num2str(median(controls_recs)) ' (' num2str(min(controls_recs)) '-' num2str(max(controls_recs)),...
    ') recordings over ' sprintf('%.0f',median(controls_range)) ' (',...
    sprintf('%.0f',min(controls_range)),...
    '-' sprintf('%.0f',max(controls_range)) ') days. ',...
    'All but ' num2str(sum(controls_missed>0)) ' subject(s) submitted at least',...
    ' one recording per week while in the study; in total only ',...
    sprintf('%.1f',100*sum(controls_alldays>=7)/length(controls_alldays)) '% of recordings',...
    ' were sent 7 or more days after the previous recording'];



%Calculate adherence
%Time on study - 24 weeks except SC05, started 11/23/20, stopped 4/7/21
%(135 days, or 19 weeks)
adherenceP = patient_recs./(24*2);
adherenceP(2) = patient_recs(2)/(19*2);
median_adP = median(adherenceP);
[rng_adP(1) rng_adP(2)] = bounds(adherenceP);
adherenceC = controls_recs./(4*2);
median_adC = mean(adherenceC);
[rng_adC(1) rng_adC(2)] = bounds(adherenceC);

text.adherence = ['Patients submitted ' sprintf('%.0f',median_adP*100) '% (range '...
    sprintf('%.0f',rng_adP(1)*100) '% to ' sprintf('%.0f',rng_adP(2)*100) '%) of'...
    ' the expected home recordings, and controls'...
    ' submitted ' sprintf('%.0f',median_adC*100) '% (range '...
    sprintf('%.0f',rng_adC(1)*100) '% to ' sprintf('%.0f',rng_adC(2)*100) '%).'];


%% Load results of Lab-evaluated 10 meter walk task
flist = dir(m2_path);
clear GA_IMU

%Load M2 parameters from the associated folder
ki = 1;
comp_tbl = [];
for i = 1:length(GA)
    clear s_path
    for j = {'tmwR','tmwL'}
        tmp2 = [GA(i).code j{:}];

        s_ind = strcmp({flist.name},tmp2);
        if sum(s_ind)==0
            disp([m2_path ' does not contain a subject folder containing ' tmp2 '\n']);
            disp('Consider running the main analysis file on all lab data');
            GA_IMU(ki).code = tmp2;
            GA_IMU(ki).WS = [];
            return
        else
            snm = flist(s_ind).name;
            s_path = [m2_path '\' snm '\' snm '_Step_Parameters'];
            s_list = dir(s_path); s_list = {s_list.name};

            step_date = regexp(s_list,...
                ['^M2parameters_imu_' snm '-(\d*)'],'tokens'); % changed M2parameters name and _ to - before the first d and removed the last 2 d's
            date_ind = find(~cellfun(@isempty,step_date));
            for n = date_ind
                if datenum(string(step_date(n)),'yyyymmdd') == datenum(GA(i).date,'m/dd/yyyy')
                    % then force date_ind to only be the date matching GA(i).date
                    date_sess = n;
                end
            end
            if isempty(date_ind)
                disp([snm ' no M2 step parameters exist for this subject']);
                GA_IMU(ki).code = tmp2;
                GA_IMU(ki).WS = [];
            else
                IMUstats = load([m2_path '\' snm '\' snm '_Step_Parameters\' s_list{date_sess}]); % changed from 1 to k ML 18-Mar-2022 because the first session was being copied twice
                IMUstats = IMUstats.stridestats;

                GA_IMU(ki).code = tmp2;
                GA_IMU(ki).date = cat(2,step_date{date_sess}{:}{1});
                GA_IMU(ki).SL = IMUstats.SL;  %/3.281; %convert to meters
                GA_IMU(ki).DU = IMUstats.DU;
                GA_IMU(ki).WS = GA_IMU(ki).SL./GA_IMU(ki).DU;
                disp(['processed IMU gait assessment data for ' tmp2 '_',...
                    GA_IMU(ki).date]);

                ki=ki+1; % moved up 3 levels
            end
        end
    end
end

% store LAB data
r = 1;
l = 1;
LAB_IMUR = cell(length(GA_IMU)/2,1);
LAB_IMUL = cell(length(GA_IMU)/2,1);
LAB_SW = cell(length(GA_IMU)/2,1);
LAB_CODE = cell(length(GA_IMU)/2,1);
for i = 1:length(GA_IMU)
    if GA_IMU(i).code(end) == 'R'
        LAB_IMUR{r} = {GA_IMU(i).WS};
        LAB_CODE(r) = regexp(GA_IMU(i).code,'S?C\d{2}','match');
        ind = strcmp({GA.code},LAB_CODE{r});
        LAB_SW{r} = GA(ind).tmw_speed; %store stopwatch speed as well
        LAB_DATE{r} = GA(ind).date;
        r = r + 1;
    elseif GA_IMU(i).code(end) == 'L'
        LAB_IMUL{l}  ={GA_IMU(i).WS};
        l = l + 1;
    end
end

pind = cellfun(@(x) isempty(x),regexp(LAB_CODE,'^C'));
cind = ~cellfun(@(x) isempty(x),regexp(LAB_CODE,'^C'));

LAB_SWp = {LAB_SW{pind}};
LAB_SWc = {LAB_SW{cind}};
LAB_IMURp = {LAB_IMUR{pind}};
LAB_IMURc = {LAB_IMUR{cind}};
LAB_IMULp = {LAB_IMUL{pind}};
LAB_IMULc = {LAB_IMUL{cind}};


%% Vicon data
VICON = cell(length(LAB_CODE),1);

% C01
C01data_1 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C01_10MW_02.csv.mat');
C01data_2 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C01_10MW_03.csv.mat');
C01data_3 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C01_10MW_04.csv.mat');
VICON{strcmp(LAB_CODE,'C01')} = [C01data_1.stride_parameters.WS',C01data_2.stride_parameters.WS',C01data_3.stride_parameters.WS'];

% C02
C02data_1 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C02_10MW_01.csv.mat');
C02data_2 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C02_10MW_02.csv.mat');
C02data_3 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C02_10MW_03.csv.mat');
VICON{strcmp(LAB_CODE,'C02')} = [C02data_1.stride_parameters.WS',C02data_2.stride_parameters.WS',C02data_3.stride_parameters.WS'];

% C03
C03data_1 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C03_10MW_03.csv.mat');
C03data_2 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C03_10MW_04.csv.mat');
C03data_3 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C03_10MW_05.csv.mat');
VICON{strcmp(LAB_CODE,'C03')} = [C03data_1.stride_parameters.WS',C03data_2.stride_parameters.WS',C03data_3.stride_parameters.WS'];

% C04
C04data_1 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C04_10MW_02.csv.mat');
C04data_2 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C04_10MW_03.csv.mat');
C04data_3 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C04_10MW_04.csv.mat');
VICON{strcmp(LAB_CODE,'C04')} = [C04data_1.stride_parameters.WS',C04data_2.stride_parameters.WS',C04data_3.stride_parameters.WS'];

% C05
C05data_1 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C05_10MW_02.csv.mat');
C05data_2 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C05_10MW_03.csv.mat');
C05data_3 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C05_10MW_04.csv.mat');
VICON{strcmp(LAB_CODE,'C05')} = [C05data_1.stride_parameters.WS',C05data_2.stride_parameters.WS',C05data_3.stride_parameters.WS'];

% C06
C06data_1 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C06_10MW_01.csv.mat');
C06data_2 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C06_10MW_02.csv.mat');
C06data_3 = load('..\LabClinicRecordings\Vicon data\StrideParameters_C06_10MW_03.csv.mat');
VICON{strcmp(LAB_CODE,'C06')} = [C06data_1.stride_parameters.WS',C06data_2.stride_parameters.WS',C06data_3.stride_parameters.WS'];

VICON = VICON((length(RCp)+1):end)';

%% Home recordings of walking speed
% Also isolate the home recording closest in date to the clinic/lab visit
clear HOME_WS HOME_WS_date ClosestHome_WS
for i = 1:length(D.u_name)
    file = dir([m2_path '\',...
        D.u_name{i} '\' D.u_name{i} '_Step_Parameters\',...
        'M2parameters*']);
    %find recording closest in time to lab date
    filenames = {file.name};
    h_date_tmp = cellfun(@(x) datenum(x,'yyyymmdd'),regexp(filenames,'(?<=-)\d{8}(?=-)','match'));
    tmp = regexp(LAB_CODE,['^' D.u_name{i}]);
    l_date_tmp = LAB_DATE(find(cellfun(@(x) ~isempty(x),tmp)));
    l_date_tmp = datenum(l_date_tmp{:});
    tmpind = find(h_date_tmp-l_date_tmp>=0);
    if isempty(tmpind)%No home recordings occurred after lab visit, find closest home recording
        [~,ind] = min(abs(h_date_tmp-l_date_tmp));
    else %Find the first recording after the lab visit
        [~,ind] = min(h_date_tmp(tmpind)-l_date_tmp);
        ind = tmpind(ind);
    end
    for j = 1:length(filenames)
        load([file(j).folder '\' file(j).name]);

        tmpdate = regexp(file(j).name,'(?<=-)\d{8}(?=-)','match');
        tmpdate = datenum(tmpdate{:},'yyyymmdd');

        HOME_WS_date{i,j} = tmpdate;
        HOME_WS{i,j} = (stridestats.SL./stridestats.DU);
        HOME_WS{i,j}(HOME_WS{i,j}<.1) = [];
        if j==ind
            ClosestHome_WS.date{i} = datestr(tmpdate);
            ClosestHome_WS.val{i} = HOME_WS{i,j};
            ClosestHome_WS.code{i} = D.u_name{i};
        end
    end

    %Sort dates
    [tmpdate sortind] = sort(cat(2,HOME_WS_date{i,:}));
    HOME_WS_date(i,1:length(sortind)) = HOME_WS_date(i,sortind);
    HOME_WS(i,1:length(sortind)) = HOME_WS(i,sortind);

end

% split into patients and controls
pind = cellfun(@(x) isempty(x),regexp(D.u_name,'^C'));
cind = ~cellfun(@(x) isempty(x),regexp(D.u_name,'^C'));

HOME_WS_date = [HOME_WS_date(pind,:); HOME_WS_date(cind,:)];
HOME_WS = [HOME_WS(pind,:); HOME_WS(cind,:)];
ClosestHome_WS = structfun(@(x) [x(pind) x(cind)],ClosestHome_WS,'UniformOutput',false);
ClosestHome_WSmean = cellfun(@(x) mean(x), ClosestHome_WS.val);


%% no significant differences between patient WS lab measures, so can combine L and R IMU
% all sessions
IMU_allp = {};
for i = 1:length(LAB_IMURp)
    IMU_allp{i} = [LAB_IMULp{i}{:},LAB_IMURp{i}{:}];
end

% results - no sigificant differences between left and right controls WS
IMU_allc = {};
for i = 1:length(LAB_IMURc)
    IMU_allc{i} = [LAB_IMULc{i}{:},LAB_IMURc{i}{:}];
end

% combine patients and controls
IMU_all = [IMU_allp,IMU_allc]; % all sessions

% find mean from each day (median was very similar)
IMU_allmean = cellfun(@(x) mean(x), IMU_all);


%% Controls lab comparison of walking speed via K-W test
% For significance bars on figure 1

h = zeros(length(IMU_allc),2);
p = zeros(length(IMU_allc),2);
for i = 1:length(LAB_SWc)
    [h(i,1),p(i,1)] = lillietest([LAB_IMURc{i}{:}'; LAB_IMULc{i}{:}']);
    [h(i,2),p(i,2)] = lillietest(VICON{i}');
end

% find mean values
LAB_SWcmean = zeros(length(LAB_SWc),1);
LAB_IMURcmean = zeros(length(LAB_SWc),1);
LAB_IMULcmean = zeros(length(LAB_SWc),1);
VICON_mean = zeros(length(LAB_SWc),1);
for i = 1:length(LAB_SWc)
    LAB_SWcmean(i) = mean(LAB_SWc{i});
    LAB_IMURcmean(i) = cellfun(@(x) mean(x),LAB_IMURc{i});
    LAB_IMULcmean(i) = cellfun(@(x) mean(x),LAB_IMULc{i});
    VICON_mean(i) = mean(VICON{i});
end

comp_tbl_mean= [];
clear PairComp
for i = 1:length(LAB_SWc)
    [pval(i),~,stats]=kruskalwallis([LAB_SWc{i}(:); LAB_IMURc{i}{:}'; LAB_IMULc{i}{:}'; VICON{i}'],...
        [ones(length(LAB_SWc{i}(:)),1);2*ones(length(LAB_IMURc{i}{:}),1);...
        2*ones(length(LAB_IMULc{i}{:}),1); 3*ones(length(VICON{i}),1)],'off');
    PairComp{i} = multcompare(stats,'display','off');
    PairComp{i}(:,1:2) = cellfun(@str2num,stats.gnames(PairComp{i}(:,1:2)));
    PairComp{i}(PairComp{i}(:,6)>.05,:) = [];
    comp_tbl_mean = [comp_tbl_mean; [mean(LAB_SWc{i}) mean([LAB_IMURc{i}{:}'; LAB_IMULc{i}{:}']),...
        mean(VICON{i}),...
        length(LAB_SWc{i}) length(LAB_IMURc{i}{:})+...
        length(LAB_IMULc{i}{:}) length(VICON{i})]];
end

ControlsLabKW = table({RCc.study_code}',comp_tbl_mean(:,1:3),comp_tbl_mean(:,4:6), pval',...
    'VariableNames',{'Study Code',...
    'Mean Walking Speed (S|M|V)','N estimates (S|M|V)','p-value'})

%% Patients lab comparison of walking speed via K-W test
% For significance bars on figure 2

clear PairComp
comp_tbl_mean= [];
for i = 1:length(LAB_SWp)
    [pval(i),~,stats] = kruskalwallis([LAB_SWp{i} LAB_IMURp{i}{:} LAB_IMULp{i}{:}],...
        [ones(length(LAB_SWp{i}),1);2*ones(length(LAB_IMURp{i}{:}),1);...
        2*ones(length(LAB_IMULp{i}{:}),1)],'off');
    PairComp{i} = multcompare(stats,'display','off');
    PairComp{i}(:,1:2) = cellfun(@str2num,stats.gnames(PairComp{i}(:,1:2)));
    PairComp{i}(PairComp{i}(:,6)>.05,:) = [];
    comp_tbl_mean = [comp_tbl_mean; [mean(LAB_SWp{i}) mean([LAB_IMURp{i}{:} LAB_IMULp{i}{:}]),...
        length(LAB_SWp{i}) length(LAB_IMURp{i}{:})+...
        length(LAB_IMULp{i}{:})]];
end

PatientsLabKW = table({RCp.study_code}',comp_tbl_mean(:,1:2),comp_tbl_mean(:,3:4), pval',...
    'VariableNames',{'Study Code',...
    'Mean Walking Speed (S|M)','N estimates (S|M)','p-value'})

%% Wilcoxon ranksum test for Home vs. Clinic IMU walking speeds

% test IMU_all and ClosestHome_WS for normality
h = zeros(length(IMU_all),2);
p = zeros(length(IMU_all),2);

for i = 1:length(IMU_all)
    [h(i,1),p(i,1)] = lillietest(IMU_all{i});
    [h(i,2),p(i,2)] = lillietest(ClosestHome_WS.val{i});
end

comp_tbl_mean= [];
clear pval pval2
for i = 1:length(IMU_all)
    %We get qualitatively similar results using ttest2
    %[~,pval(i),~,stats] = ttest2(IMU_all{i}, ClosestHome_WS{i});
    [pval(i),~,stats] = ranksum(IMU_all{i}, ClosestHome_WS.val{i});
    comp_tbl_mean = [comp_tbl_mean; [mean(IMU_all{i}) mean(ClosestHome_WS.val{i}),...
        length(IMU_all{i}) length(ClosestHome_WS.val{i})]];
end

LabHomeRS = table({RC.study_code}',comp_tbl_mean(:,1:2),comp_tbl_mean(:,3:4), pval',...
    'VariableNames',{'Study Code',...
    'Mean Walking Speed (L|H)','N estimates (L|H)','p-value'})

%% Linear Mixed Model for home vs clinic speed
% want to compare clinic and home imu WS for patients and controls
% combine data into a table
% columns: subject, clinic/home (0/1), WS
% rows: data

code_all = [3, 5, 7, 8, 9, 10, 11, 12, 13, 21, 22, 23, 24, 25, 26];
%location identifier
H = 1;
C = 0;

% linear model based on mean values per session
lme_names_2 = [];
sub_type2 = [];
lme_wsid_2 = [];
lme_wsdata_2 = [];
for i = 1:length(IMU_allmean)
    lme_names_2 = [lme_names_2;code_all(i);code_all(i)];
    lme_wsid_2 = [lme_wsid_2;C;H];
    lme_wsdata_2 = [lme_wsdata_2;IMU_allmean(i);ClosestHome_WSmean(i)];
    if sum(i == [1:9])>0 % patient
        ST = 0;
    else
        ST = 1;
    end
    sub_type2 = [sub_type2;ST;ST];
end

% form into a table
headings_mean = {'Code','Identifier','WS','SubjectType'};
lme_tbl_data_2 = [lme_names_2,lme_wsid_2,lme_wsdata_2,sub_type2];
lme_tbl_2 = array2table(lme_tbl_data_2,'VariableNames',headings_mean);
lme_2 = fitlm(lme_tbl_2,'WS ~ Identifier + SubjectType')


%% Walking speed over time

code_pat = LAB_CODE(1:9);

ws_data_all = [];
ws_names_all = [];
ws_days_all = [];
ws_mean_all = [];
days_mean_all = [];
names_mean_all = [];
for i = 1:length(code_pat)
    tmpdate = cellfun(@(x,y) repmat(y,length(x),1),HOME_WS(i,:),HOME_WS_date(i,:),'UniformOutput',false);
    [tmpdate, sortdate] = sort(cat(1,tmpdate{:}));
    ws_days_all = [ws_days_all;tmpdate-tmpdate(1)];
    tmpws = cat(2,HOME_WS{i,:})';
    ws_data_all = [ws_data_all;tmpws(sortdate)];
    ws_names_all = [ws_names_all;repmat(code_pat(i),length(tmpws),1)];

    %Daily means
    tmpdatemean = [HOME_WS_date{i,:}];
    [tmpdatemean, sortdatemean] = sort(tmpdatemean);
    tmpwsmean = cellfun(@mean,HOME_WS(i,:))';
    tmpwsmean(isnan(tmpwsmean)) = [];
    ws_mean_all = [ws_mean_all;tmpwsmean(sortdatemean)];
    days_mean_all = [days_mean_all;tmpdatemean'-tmpdatemean(1)];
    names_mean_all = [names_mean_all;repmat(code_pat(i),length(tmpwsmean),1)];
end

% form into a table
headings = {'Code','Days','WalkingSpeed'};
ws_tbl_data = {ws_names_all,ws_days_all,ws_data_all};
ws_tbl = table(ws_names_all,ws_days_all,ws_data_all,'VariableNames',headings);

% table for means
headings_mean = {'Code','Days','WalkingSpeed'};
ws_tbl_means_data = {names_mean_all,days_mean_all,ws_mean_all};
ws_tbl_means = table(names_mean_all,days_mean_all,ws_mean_all,'VariableNames',headings_mean);

[a,b,c] = unique(ws_tbl_means_data{1});
firstws = ws_tbl_means_data{3}(b);
lme_hws = fitlm(ws_tbl_means,'WalkingSpeed ~ Days')

text.wsovertime = ['Patients'' initial home walking speeds ranged from ',...
    sprintf('%.2f',min(firstws)) '-' sprintf('%.2f',max(firstws)) 'm/s. ',...
    sprintf('%.0f',sum(firstws<.7)) ' patients started with a walking speed',...
    ' less that .7 m/s,' sprintf('%.0f',sum(firstws>=.7 & firstws<=1)),...
    ' between .7 and 1, and ' sprintf('%.0f',sum(firstws>1)) ' patients',...
    ' greater than 1 m/s.',...
    ' Change in walking speed was ' sprintf('%.4f',lme_hws.Coefficients.Estimate(2)),...
    ' m/s/day (p = ' sprintf('%.3f',lme_hws.Coefficients.pValue(2)) ').'];


%% Linear model, ws over time change compared to assist device changes
% for subject SC7, SC8, and SC11 who changed assist devices
% use code, days, and ws data for specific subjects from above cell

newlist = {'SC07','SC08','SC11'};
assistdata_ch = cell(1,length(newlist));

clear fnames fdatei folder fold_dates
tmp_assist = cell(1,length(assist_sub));
for i = 1:length(newlist)
    % take ws over time data with days and subject codes
    assistdata_ch{i} = cellfun(@(x) x(strcmp(ws_tbl_data{1},newlist{i}),:),ws_tbl_data,'UniformOutput',false);

    % add on assist data, using assist variable from earlier cell
    fldr = dir([m2_path '\' newlist{i} '\' newlist{i} '_Step_Parameters\M2*']);
    % make sure files are in order
    fnames = {fldr.name};
    fdatei = cellfun(@(x) strfind(x,'-'),fnames,'UniformOutput',false);
    clear fdate
    for n = 1:length(fdatei)
        tmp = cell2mat(fdatei(n));
        fdate(n) = datenum(fnames{n}(tmp(1)+1:tmp(2)-1),'yyyymmdd');
    end
    [~,d] = sort(fdate);
    k = 2;
    for j = 1:length(fnames)
        clear tmpdata tmpfile
        tmpfile = load([fldr(d(j)).folder '\' fldr(d(j)).name]);


        tmpdata = (tmpfile.stridestats.SL./tmpfile.stridestats.DU)';
        tmp_assist{i} = [tmp_assist{i};repmat(assist_sub{i}(k),length(tmpdata),1)];
        k = k + 1;
    end
end

% append tmp_assist to assistdata
for i = 1:length(tmp_assist)
    assistdata_ch{i}{4} = tmp_assist{i};
end

% stack subject data
assist_tbl_data = [];
for i = 1:length(assistdata_ch)
    assist_tbl_data = [assist_tbl_data;assistdata_ch{i}];
end

% form into a table
headings = {'Code','Days','WalkingSpeed','AssistDevice'};
assist_tbl = array2table(assist_tbl_data,'VariableNames',headings);

%Create table of means
assist_data_mean = [];
for i = 1:length(newlist)
    tmpws = cellfun(@mean,HOME_WS(strcmp(LAB_CODE,newlist{i}),:));
    tmpws(isnan(tmpws)) = [];
    tmpdate = cellfun(@mean,HOME_WS_date(strcmp(LAB_CODE,newlist{i}),:));
    tmpdate(isnan(tmpdate)) = [];
    tmpdate = tmpdate-tmpdate(1);
    tmpassist = assist_sub{i}(2:end)>0;
    assist_data_mean = [assist_data_mean;...
        [repmat(i,length(tmpws),1) tmpdate' tmpws' tmpassist ]];

end
assist_tbl_mean = table(assist_data_mean(:,1),assist_data_mean(:,2),...
    assist_data_mean(:,3),assist_data_mean(:,4),...
    'VariableNames',{'Subject','Days','WalkingSpeed','AssistDevice'});

clear lme_ast
for i = 1:length(newlist)
    ind = assist_data_mean(:,1)==i;
    lme_ast{i} = fitlm(assist_tbl_mean(ind,:),'WalkingSpeed ~ AssistDevice*Days');
end

%Calculate instantaneous change in walking speed following use of walking device
%SC07
%From Eqn 5
%v_none = phi0 + phi1*Days
v_none7 = lme_ast{1}.Coefficients.Estimate(1) + ...
    lme_ast{1}.Coefficients.Estimate(2)*assist_data_mean(assist_data_mean(:,1)==1 & [0; diff(assist_data_mean(:,4))>0],2);
%v_ast = phi0 + phi2 + (phi1 + phi3) * Days
v_ad7 = (lme_ast{1}.Coefficients.Estimate(1) + lme_ast{1}.Coefficients.Estimate(3)) + ...
    (lme_ast{1}.Coefficients.Estimate(2) + lme_ast{1}.Coefficients.Estimate(4))*...
    assist_data_mean(assist_data_mean(:,1)==1 & [0; diff(assist_data_mean(:,4))>0],2);
v_diff7 = v_none7-v_ad7;
%SC08
v_none8 = lme_ast{2}.Coefficients.Estimate(1) + ...
    lme_ast{2}.Coefficients.Estimate(2)*assist_data_mean(assist_data_mean(:,1)==2 & [0; diff(assist_data_mean(:,4))>0],2);
v_ad8 = (lme_ast{2}.Coefficients.Estimate(1) + lme_ast{2}.Coefficients.Estimate(3)) + ...
    (lme_ast{2}.Coefficients.Estimate(2) + lme_ast{2}.Coefficients.Estimate(4))*...
    assist_data_mean(assist_data_mean(:,1)==2 & [0; diff(assist_data_mean(:,4))>0],2);
v_diff8 = v_none8-v_ad8;
%SC11
v_none11 = lme_ast{3}.Coefficients.Estimate(1) + ...
    lme_ast{3}.Coefficients.Estimate(2)*assist_data_mean(assist_data_mean(:,1)==3 & [0; diff(assist_data_mean(:,4))>0],2);
v_ad11 = (lme_ast{3}.Coefficients.Estimate(1) + lme_ast{3}.Coefficients.Estimate(3)) + ...
    (lme_ast{3}.Coefficients.Estimate(2) + lme_ast{3}.Coefficients.Estimate(4))*...
    assist_data_mean(assist_data_mean(:,1)==3 & [0; diff(assist_data_mean(:,4))>0],2);
v_diff11 = v_none11-v_ad11;


%% Plot figures
statistics_figs_v4(LAB_SW,IMU_all,VICON,ClosestHome_WS,HOME_WS,HOME_WS_date,assist_tbl_data,lme_ast)
