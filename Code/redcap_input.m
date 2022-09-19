%This function reads the .csv files on the path specified and processes
%them for relevant demographic, clinical, and assessment data.
function [RCdata,RCind,record_id_vec] = redcap_input(redcap_path,proc_path,Vbl)
disp '    Begin processing Redcap data';

bigOrg = redcap_path.bigOrg;
rawall = redcap_path.rawall;
indiv = redcap_path.indiv;
group = redcap_path.group;

path_csv = proc_path.csv;
path_fig = proc_path.fig;
path_results = proc_path.results;
path_datafig = proc_path.datafig;
path_pathfig = proc_path.pathfig;
path_preprocfig = proc_path.preprocfig;
path_imu = proc_path.imu;
path_step = proc_path.step;
path_line = proc_path.line;


% begin file organization
if ~exist(bigOrg, 'dir')       % check to make sure file org folder exists
    mkdir(bigOrg);              % if it doesn't, create it
end

if ~exist([bigOrg,filesep,group], 'dir')
    mkdir ([bigOrg,filesep,group])
end

if ~exist([bigOrg,filesep,group,filesep,path_fig], 'dir')
    mkdir ([bigOrg,filesep,group,filesep,path_fig])
end
%%%%%%

if ~exist([bigOrg,filesep,indiv], 'dir')
    mkdir ([bigOrg,filesep,indiv])
end

if length(Vbl.pathtoredcap)==2 %Home Data
    rc_types = {'P','C'}; %Patients and controls have separate redcap files
    [~,~,rawC] = xlsread(Vbl.pathtoredcap{1});
    [~,~,rawP] = xlsread(Vbl.pathtoredcap{2});
else
    rc_types = {'A'}; %Patients and Controls in one file
    [~,~,rawA] = xlsread(Vbl.pathtoredcap{1});
end
    
%Loop through patients and controls
gp_ct = 0;
for st = rc_types
    raw = eval(['raw' st{:}]);
    
    % enrollment data:
    RCind.recordid = find(strcmp(raw(1,:),'record_id'));           % study code
    RCind.recordtype = find(strcmp(raw(1,:),'redcap_repeat_instrument'));         % type of record
    RCind.repeat = find(strcmp(raw(1,:),'redcap_repeat_instance'));             % how many repeats
    RCind.enrolldate = find(strcmp(raw(1,:),'enroll_date'));         % person's enrollment date
    RCind.age = find(strcmp(raw(1,:),'age'));                % person's age
    RCind.sex = find(strcmp(raw(1,:),'sex'));                % person's sex
    RCind.phone = find(strcmp(raw(1,:),'phone_type'));              % person's phone type
    RCind.enrollment = find(strcmp(raw(1,:),'enrollment_complete'));         % incomplete = 0, complete = 1
    
    % historical data:
    RCind.h_visitdate = find(strcmp(raw(1,:),'hist_visit_date'));         % date of last clinic visit
    RCind.h_height = find(strcmp(raw(1,:),'hist_height'));
    RCind.h_weight = find(strcmp(raw(1,:),'hist_weight'));
    RCind.h_falls = find(strcmp(raw(1,:),'hist_reported_falls'));
    RCind.h_frsdate = find(strcmp(raw(1,:),'hist_alsfrs_date'));
    RCind.h_frs = find(strcmp(raw(1,:),'hist_frs_1')):find(strcmp(raw(1,:),'hist_frs_12'));
    
    RCind.h_strdate = find(strcmp(raw(1,:),'hist_strength_date'));           % date of last strength test
    RCind.h_strhip_r = find(strcmp(raw(1,:),'hist_str_hip_flex_r'));          % last strength test: right hip flex
    RCind.h_strkflex_r = find(strcmp(raw(1,:),'hist_str_knee_flex_r'));        % last strength test: right knee flex
    RCind.h_strkext_r = find(strcmp(raw(1,:),'hist_str_knee_ext_r'));         % last strength test: right knee extension
    RCind.h_strfoot_r = find(strcmp(raw(1,:),'hist_str_foot_dorsiflex_r'));         % last strength test: right foot dorsiflex
    RCind.h_strhip_l = find(strcmp(raw(1,:),'hist_str_hip_flex_l'));          % last strength test: left hip flex
    RCind.h_strkflex_l = find(strcmp(raw(1,:),'hist_str_knee_flex_l'));        % last strength test: left knee flex
    RCind.h_strkext_l = find(strcmp(raw(1,:),'hist_str_knee_ext_l'));         % last strength test: left knee extension
    RCind.h_strfoot_l = find(strcmp(raw(1,:),'hist_str_foot_dorsiflex_l'));         % last strength test: left foot dorsiflex
    
    RCind.h_umndate = find(strcmp(raw(1,:),'hist_umn_date'));           % date of last reflex test
    RCind.h_umnpat_r = find(strcmp(raw(1,:),'hist_umn_pat_r'));          % last reflex test: right patellar
    RCind.h_umnach_r = find(strcmp(raw(1,:),'hist_umn_ach_r'));        % last reflex test: right achilles
    RCind.h_umnpat_l = find(strcmp(raw(1,:),'hist_umn_pat_l'));         % last reflex test: left patellar
    RCind.h_umnach_l = find(strcmp(raw(1,:),'hist_umn_ach_l'));         % last reflex test: left achilles
    
    
    RCind.h_fvcdate = find(strcmp(raw(1,:),'hist_fvc_date'));           % date of last Forced Vital Capacity test
    RCind.h_fvcliters = find(strcmp(raw(1,:),'hist_fvc_liters'));         % last FVC liters of air
    RCind.h_fvcpp = find(strcmp(raw(1,:),'hist_fvc_pp'));             % last FVC percent lung capacity
    RCind.h_note = find(strcmp(raw(1,:),'hist_note'));              % historical data note
    RCind.h_complete = find(strcmp(raw(1,:),'historical_clinical_data_complete'));          % complete history data
    
    % clinical data:
    RCind.visitdate = find(strcmp(raw(1,:),'visit_date'));          % date of clinic visit
    RCind.height = find(strcmp(raw(1,:),'height'));             % patient's height
    RCind.weight = find(strcmp(raw(1,:),'weight'));            % patient's weight
    RCind.falls = find(strcmp(raw(1,:),'reported_falls'));               % patient's past reported falls
    RCind.frsdate = find(strcmp(raw(1,:),'alsfrs_date'));           % date of ALSFRS test
    RCind.frs = find(strcmp(raw(1,:),'frs_1')):find(strcmp(raw(1,:),'frs_12'));            % ALSFRS data
    RCind.frspeg = find(strcmp(raw(1,:),'frs_peg'));            %
    
    RCind.strdate = find(strcmp(raw(1,:),'strength_date'));           % date of strength test
    RCind.strhip_r = find(strcmp(raw(1,:),'strength_hip_flex_r'));          % strength test: right hip flex
    RCind.strkflex_r = find(strcmp(raw(1,:),'strength_knee_flex_r'));        % strength test: right knee flex
    RCind.strkext_r = find(strcmp(raw(1,:),'strength_knee_ext_r'));         % strength test: right knee extension
    RCind.strfoot_r = find(strcmp(raw(1,:),'strength_foot_dorsiflex_r'));         % strength test: right foot dorsiflex
    RCind.strhip_l = find(strcmp(raw(1,:),'strength_hip_flex_l'));          % strength test: left hip flex
    RCind.strkflex_l = find(strcmp(raw(1,:),'strength_knee_flex_l'));        % strength test: left knee flex
    RCind.strkext_l = find(strcmp(raw(1,:),'strength_knee_ext_l'));         % strength test: left knee extension
    RCind.strfoot_l = find(strcmp(raw(1,:),'strength_foot_dorsiflex_l'));         % strength test: left foot dorsiflex
    
    RCind.umndate = find(strcmp(raw(1,:),'umn_date'));           % date of reflex test
    RCind.umnpat_r = find(strcmp(raw(1,:),'umn_pat_r'));          %  reflex test: right patellar
    RCind.umnach_r = find(strcmp(raw(1,:),'umn_ach_r'));        %  reflex test: right achilles
    RCind.umnpat_l = find(strcmp(raw(1,:),'umn_pat_l'));         %  reflex test: left patellar
    RCind.umnach_l = find(strcmp(raw(1,:),'umn_ach_l'));         %  reflex test: left achilles
    
    RCind.fvcdate = find(strcmp(raw(1,:),'fvc_date'));           % date of Forced Vital Capacity test
    RCind.fvcliters = find(strcmp(raw(1,:),'fvc_liters'));         % FVC liters of air
    RCind.fvcpp = find(strcmp(raw(1,:),'fvc_pp'));             % FVC percent lung capacity
    
    RCind.meddate = find(strcmp(raw(1,:),'medication_date'));           % date of filling out list of medications
    RCind.medname = find(strcmp(raw(1,:),'medication'));           % medication name
    
    RCind.fesdate = find(strcmp(raw(1,:),'fes_date'));           % date of Falls Efficacy Scale
    RCind.fes = find(strcmp(raw(1,:),'fes_clean')):find(strcmp(raw(1,:),'fes_social'));            % FES values
    RCind.visitnote = find(strcmp(raw(1,:),'visit_note'));         % visit note
    RCind.clinicdata = find(strcmp(raw(1,:),'clinical_data_complete'));        % incomplete clinic data = 0, complete = 1
    
    % gait and motor assessment:
    RCind.gadate = find(strcmp(raw(1,:),'ga_date'));            % date of gait and motor assessment
    RCind.leglengthl = find(strcmp(raw(1,:),'leg_length_left'));        % person's LEFT leg length
    RCind.leglengthr = find(strcmp(raw(1,:),'leg_length_right'));        % person's RIGHT leg length
    RCind.mobdev = find(strcmp(raw(1,:),'mobility_devices'));            % mobility devices used
    RCind.brace = find(strcmp(raw(1,:),'braces_used'));             % braces used
    RCind.sensorr = find(strcmp(raw(1,:),'sensorid_rf'));           % RIGHT foot's sensor ID
    RCind.sensorl = find(strcmp(raw(1,:),'sensorid_lf'));           % LEFT foot's sensor ID
    RCind.sensort = find(strcmp(raw(1,:),'sensorid_t'));           % TORSO's sensor ID
    
    RCind.walkdata = find(strcmp(raw(1,:),'walk_data'));          % .zip file of 10-m walk data
    RCind.walk = find(strcmp(raw(1,:),'walk_1_t')):find(strcmp(raw(1,:),'walk_3_t'));           % times of trials 1-3, respectively
    RCind.tugdata = find(strcmp(raw(1,:),'tug_data'));           % .zip file of TUG data
    RCind.tug = find(strcmp(raw(1,:),'tug_1_t')):find(strcmp(raw(1,:),'tug_3_t'));            % times of trials 1-3, respectively
    
    RCind.berg = find(strcmp(raw(1,:),'berg_1')):find(strcmp(raw(1,:),'berg_14'));           % Berg balance scale criteria
    RCind.bergdata = find(strcmp(raw(1,:),'berg_data'));          % .zip file of sensor data from Berg test
    RCind.gamnotes = find(strcmp(raw(1,:),'gait_assessment_notes'));          % notes on gait and motor assessment
    RCind.gamcomplete = find(strcmp(raw(1,:),'gait_and_motor_function_assessment_complete'));       % gait and motor assessment incomplete = 0
    
    RCind.studycomplete = find(strcmp(raw(1,:),'complete_study'));      % study completed
    RCind.withdraw = find(strcmp(raw(1,:),'withdraw_reason'));           % reason for withdrawing from the study
    RCind.comments = find(strcmp(raw(1,:),'study_comments'));           % study comments
    RCind.completioncomplete = find(strcmp(raw(1,:),'completion_data_complete'));
    
    clear record_id_vec
    for i=2:size(raw,1)
        record_id = raw{i,RCind.recordid};     % indidivual's study code
        record_id_vec{i-1} = {num2str(record_id)};        % initialize outside loop
        record_type = raw{i,RCind.recordtype}; % can be 'enrollment','clinical_data',or 'gait and motor function assessment'
        record_repeat = raw{i,RCind.repeat};   % the number of iterations for one person (including current one)
        if isnan(record_type)
            record_type = 'enrollment';
            record_repeat = 1;
        end
        
        ptcount = gp_ct+size(unique(cellfun(@cell2mat,record_id_vec, 'UniformOutput',false)),2);  % running total # of pt's
        
        switch record_type
            case 'enrollment'
                tmp = raw{i,RCind.recordid};
                if isnumeric(tmp)
                    tmp = sprintf('SC%02d',tmp);
                end
                RCdata(ptcount).study_code = tmp;
                RCdata(ptcount).enroll_date = raw{i,RCind.enrolldate};
                RCdata(ptcount).e_age = raw{i,RCind.age};
                if strcmpi(raw(i,RCind.sex),'female')
                    raw{i,RCind.sex} = 0;
                elseif strcmpi(raw(i,RCind.sex),'male')
                    raw{i,RCind.sex} = 1;
                else
                    raw{i,RCind.sex} = NaN;
                end
                RCdata(ptcount).e_sex{record_repeat} = raw{i,RCind.sex};
                RC(ptcount).e_phone = raw{i,RCind.phone};
                % complete enrollment
                RCdata(ptcount).e_complete(record_repeat) = raw(i,RCind.enrollment);
            case 'historical_clinical_data'
                % Historical Visit:
                
                RCdata(ptcount).h_visit_date(record_repeat) = raw(i,RCind.h_visitdate);
                RCdata(ptcount).h_height(record_repeat) = raw(i,RCind.h_height);
                RCdata(ptcount).h_weight(record_repeat) = raw(i,RCind.h_weight);
                % ALSFRS test
                RCdata(ptcount).h_alsfrs_date(record_repeat) = raw(i,RCind.h_frsdate);
                RCdata(ptcount).h_alsfrs(record_repeat) = {raw(i,RCind.h_frs)};
                % Strength test right limbs
                RCdata(ptcount).h_str_date(record_repeat) = raw(i,RCind.h_strdate);
                RCdata(ptcount).h_str_hip_r(record_repeat) = raw(i,RCind.h_strhip_r);
                RCdata(ptcount).h_str_k_flex_r(record_repeat) = raw(i,RCind.h_strkflex_r);
                RCdata(ptcount).h_str_k_ext_r(record_repeat) = raw(i,RCind.h_strkext_r);
                RCdata(ptcount).h_str_foot_r(record_repeat) = raw(i,RCind.h_strfoot_r);
                % Strength test left limbs
                RCdata(ptcount).h_str_hip_l(record_repeat) = raw(i,RCind.h_strhip_l);
                RCdata(ptcount).h_str_k_flex_l(record_repeat) = raw(i,RCind.h_strkflex_l);
                RCdata(ptcount).h_str_k_ext_l(record_repeat) = raw(i,RCind.h_strkext_l);
                RCdata(ptcount).h_str_foot_l(record_repeat) = raw(i,RCind.h_strfoot_l);
                % Upper motor neuron test right limbs
                RCdata(ptcount).h_umn_date(record_repeat) = raw(i,RCind.h_umndate);
                RCdata(ptcount).h_umn_pat_r(record_repeat) = raw(i,RCind.h_umnpat_r);
                RCdata(ptcount).h_umn_ach_r(record_repeat) = raw(i,RCind.h_umnach_r);
                % Upper motor neuron test left limbs
                RCdata(ptcount).h_umn_pat_l(record_repeat) = raw(i,RCind.h_umnpat_l);
                RCdata(ptcount).h_umn_ach_l(record_repeat) = raw(i,RCind.h_umnach_l);
                % Forced vital capacity test
                RCdata(ptcount).h_fvc_date(record_repeat) = raw(i,RCind.h_fvcdate);
                RCdata(ptcount).h_fvc_liters(record_repeat) = raw(i,RCind.h_fvcliters);
                RCdata(ptcount).h_fvc_pp(record_repeat) = raw(i,RCind.h_fvcpp);
                % complete clinical data
                RCdata(ptcount).h_complete(record_repeat) = {raw(i,RCind.h_complete)};
                
            case 'clinical_data' % added this case because it did not exist before
                % Current Visit:
                
                if ~isnan(cell2mat(raw(i,RCind.visitdate))) %if date empty, skip over entries
                    RCdata(ptcount).c_visit_date(record_repeat) = raw(i,RCind.visitdate);
                    RCdata(ptcount).c_height(record_repeat) = raw(i,RCind.height);
                    RCdata(ptcount).c_weight(record_repeat) = raw(i,RCind.weight);
                end
                
                % ALSFRS test
                if ~isnan(cell2mat(raw(i,RCind.frsdate))) %if date empty, skip over entries
                    RCdata(ptcount).c_alsfrs_date(record_repeat) = raw(i,RCind.frsdate);
                    RCdata(ptcount).c_alsfrs(record_repeat) = {raw(i,RCind.frs)};
                    RCdata(ptcount).c_frs_peg(record_repeat) = raw(i,RCind.frspeg);
                end
                
                % Strength test right limbs
                if ~isnan(cell2mat(raw(i,RCind.frsdate))) %if date empty, skip over entries
                    RCdata(ptcount).c_str_date(record_repeat) = raw(i,RCind.strdate);
                    RCdata(ptcount).c_str_hip_r(record_repeat) = raw(i,RCind.strhip_r);
                    RCdata(ptcount).c_str_k_flex_r(record_repeat) = raw(i,RCind.strkflex_r);
                    RCdata(ptcount).c_str_k_ext_r(record_repeat) = raw(i,RCind.strkext_r);
                    RCdata(ptcount).c_str_foot_r(record_repeat) = raw(i,RCind.strfoot_r);
                    % Strength test left limbs
                    RCdata(ptcount).c_str_hip_l(record_repeat) = raw(i,RCind.strhip_l);
                    RCdata(ptcount).c_str_k_flex_l(record_repeat) = raw(i,RCind.strkflex_l);
                    RCdata(ptcount).c_str_k_ext_l(record_repeat) = raw(i,RCind.strkext_l);
                    RCdata(ptcount).c_str_foot_l(record_repeat) = raw(i,RCind.strfoot_l);
                end
                
                % Upper motor neuron test right limbs
                if ~isnan(cell2mat(raw(i,RCind.umndate))) %if date empty, skip over entries
                    RCdata(ptcount).c_umn_date(record_repeat) = raw(i,RCind.umndate);
                    RCdata(ptcount).c_umn_pat_r(record_repeat) = raw(i,RCind.umnpat_r);
                    RCdata(ptcount).c_umn_ach_r(record_repeat) = raw(i,RCind.umnach_r);
                    % Upper motor neuron test left limbs
                    RCdata(ptcount).c_umn_pat_l(record_repeat) = raw(i,RCind.umnpat_l);
                    RCdata(ptcount).c_umn_ach_l(record_repeat) = raw(i,RCind.umnach_l);
                end
                
                % Forced vital capacity test
                if ~isnan(cell2mat(raw(i,RCind.fvcdate))) %if date empty, skip over entries
                    RCdata(ptcount).c_fvc_date(record_repeat) = raw(i,RCind.fvcdate);
                    RCdata(ptcount).c_fvc_liters(record_repeat) = raw(i,RCind.fvcliters);
                    RCdata(ptcount).c_fvc_pp(record_repeat) = raw(i,RCind.fvcpp);
                end
                
                % Data on medications
                if ~isnan(cell2mat(raw(i,RCind.meddate))) %if date empty, skip over entries
                    RCdata(ptcount).c_med_date(record_repeat) = raw(i,RCind.meddate);
                    RCdata(ptcount).c_med_name(record_repeat) = raw(i,RCind.medname);
                end
                
                % Falls efficacy scale
                if ~isnan(cell2mat(raw(i,RCind.fesdate))) %if date empty, skip over entries
                    RCdata(ptcount).c_fes_date(record_repeat) = raw(i,RCind.fesdate);
                    RCdata(ptcount).c_fes(record_repeat) = {raw(i,RCind.fes)};
                end
                
                % complete clinical data
                RCdata(ptcount).c_complete(record_repeat) = {raw(i,RCind.clinicdata)};
            case 'gait_and_motor_function_assessment'
                if ~isnan(cell2mat(raw(i,RCind.gadate))) %if date empty, skip over entries
                    RCdata(ptcount).gam_ga_date(record_repeat) = raw(i,RCind.gadate);
                    RCdata(ptcount).gam_leg_length_L(record_repeat) = raw(i,RCind.leglengthl);
                    RCdata(ptcount).gam_leg_length_R(record_repeat) = raw(i,RCind.leglengthr);
                    RCdata(ptcount).gam_mob_dev(record_repeat) = raw(i,RCind.mobdev);
                    RCdata(ptcount).gam_brace(record_repeat) = raw(i,RCind.brace);
                    RCdata(ptcount).gam_sensor_R(record_repeat) = raw(i,RCind.sensorr);
                    RCdata(ptcount).gam_sensor_L(record_repeat) = raw(i,RCind.sensorl);
                    RCdata(ptcount).gam_sensor_T(record_repeat) = raw(i,RCind.sensort);
                    % 10-m walk and TUG data
                    RCdata(ptcount).gam_walk_data(record_repeat) = raw(i,RCind.walkdata);
                    RCdata(ptcount).gam_walk(record_repeat) = {raw(i,RCind.walk)};
                    RCdata(ptcount).gam_tug_data(record_repeat) = raw(i,RCind.tugdata);
                    RCdata(ptcount).gam_tug(record_repeat) = {raw(i,RCind.tug)};
                    % Berg balance scale
                    RCdata(ptcount).gam_berg(record_repeat) = {raw(i,RCind.berg)};
                    RCdata(ptcount).gam_berg_data(record_repeat) = raw(i,RCind.bergdata);
                    % Notes and completion of Gait and Motor Assessment
                    RCdata(ptcount).gam_notes(record_repeat) = raw(i,RCind.gamnotes);
                    RCdata(ptcount).gam_complete(record_repeat) = raw(i,RCind.gamcomplete);
                end
        end
    end
    gp_ct = ptcount;
    record_id_vec = unique(cellfun(@cell2mat,record_id_vec,'UniformOutput',false));
    RCdata(all(cell2mat(arrayfun(@(x) structfun(@isempty, x ),RCdata,'UniformOutput',false)),1)) = [];
end
disp 'Redcap data up to date';
% end redcap inputs


%Create folder structure for saving inidividual results
totalSubjects = length(RCdata); % returns total # of patients from redcap data
for i=1:totalSubjects
    record_id = RCdata(i).study_code;     % indidivual's study code
    
    indivFile = [bigOrg,filesep,indiv,filesep,record_id];
    
    if ~exist(indivFile, 'dir')
        mkdir (indivFile)
    end

    if ~exist([indivFile,filesep,record_id,path_fig], 'dir')
        mkdir([indivFile,filesep,record_id,path_fig]);
    end
    
    if ~exist([indivFile,filesep,record_id,path_fig,filesep,record_id,path_datafig], 'dir')
        mkdir([indivFile,filesep,record_id,path_fig,filesep,record_id,path_datafig]);
    end
    
    if ~exist([indivFile,filesep,record_id,path_fig,filesep,record_id,path_pathfig], 'dir')
        mkdir([indivFile,filesep,record_id,path_fig,filesep,record_id,path_pathfig]);
    end
    
    if ~exist([indivFile,filesep,record_id,path_fig,filesep,record_id,path_preprocfig], 'dir')
        mkdir([indivFile,filesep,record_id,path_fig,filesep,record_id,path_preprocfig]);
    end
    
    if ~exist([indivFile,filesep,record_id,path_imu], 'dir')
        mkdir([indivFile,filesep,record_id,path_imu]);
    end
    
    if ~exist([indivFile,filesep,record_id,path_step], 'dir')
        mkdir([indivFile,filesep,record_id,path_step]);
    end
    
    if ~exist([indivFile,filesep,record_id,path_line], 'dir')
        mkdir([indivFile,filesep,record_id,path_line]);
    end
end
disp 'Subject folders up to date';
end
