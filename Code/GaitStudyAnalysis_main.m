%% GaitStudyAnalysis_main.m
%This is the main file used for analysis of data, and the production of
%results and figures in "Longitudinal analysis of spatial-temporal gait 
%parameters for ALS patients using remotely collected foot-worn IMU data"  
%
%Authors contributing to this code: Meghan Lukac, Hannah Luben, 
%Anne E. Martin, and Andrew Geronimo 
%
%For questions regarding this code, use the correspondence address provided
%in the publication.
%
%This code was run and tested with the following package versions
%     {'MATLAB'                                 }    {'9.11'}
%     {'Signal Processing Toolbox'              }    {'8.7' }
%     {'Statistics and Machine Learning Toolbox'}    {'12.2'}
%     {'Sensor Fusion and Tracking Toolbox'     }    {'2.2' }
%     {'Navigation Toolbox'                     }    {'2.1' }


clear all
close all

%% Set up folder structure
redcap_path.bigOrg = 'Processed_Data';
redcap_path.rawall = 'Raw_Data';
redcap_path.group = 'Group_Comparisons';
redcap_path.indiv = 'Individual_Patient_Data';

proc_path.csv = '_CSV_Data';
proc_path.fig = '_Figures';
proc_path.results = '_Results';
proc_path.datafig = '_Data_Figs';
proc_path.pathfig = '_Path_Figs';
proc_path.preprocfig = '_Pre_Proc_Figs';
proc_path.imu = '_Processed_IMU';
proc_path.step = '_Step_Parameters';
proc_path.line = '_Straight_Line_Segs';

Vbl.accHeader = 'ACC';                
Vbl.gyrHeader = 'GYR';                
Vbl.figuresOn = 0;
Vbl.g_conversion = 9.81;

%% Prompt to process home or lab recordings
%Point to csv files and redcap data associated with chosen recordings
rec_i = input('Process home (0) or lab (1) recordings?\n');

if rec_i==0
    Vbl.pathtodata = ['..' filesep 'HomeRecordings' filesep];
    Vbl.pathtoredcap = {['..' filesep 'Redcap_data' filesep 'redcap_ControlData.csv'];
        ['..' filesep 'Redcap_data' filesep 'redcap_PatientData.csv']};
elseif rec_i==1
    Vbl.pathtodata = ['..' filesep 'LabClinicRecordings' filesep];
    Vbl.pathtoredcap = {['..' filesep 'Redcap_data' filesep 'redcap_LabData.csv']};
else
    disp('Invalid selection');
    return;
end

%% Load REDCap data
[RCdata,RCind,record_id_vec] = redcap_input(redcap_path,proc_path,Vbl);

%% Prompt for which data to process
F = dir([Vbl.pathtodata '*.csv']);
Fnames = {F.name};
Snames = regexpi(Fnames,'s?c\d+(tmw[LR])?','match');
Snames = unique(cellfun(@upper,[Snames{:}],'UniformOutput',false));

for i = 1:length(Snames)
    fprintf(['(' num2str(i) ') ' Snames{i} '\n'])
end
subj_i = input('Select the number of the subject to analyze.  Type 0 for all. \n');
if subj_i==0
   Vbl.subj = Snames;
elseif subj_i<=length(Snames)
    Vbl.subj = Snames(subj_i);
else
    disp('Invalid selection');
    return;
end

%% Process IMU excel files
%Results in subject-specific DATA.mat files saved in the "Home" or "LabClinic"
%folders.  Existing DATA files will not be overwritten.  To rerun
%process_imu_data on your selected set, delete the associated DATA files.
%This can be done manually or using clean_results.m
[numericdata, metadata] = process_imu_data(redcap_path,proc_path,Vbl);

%% Manually choose walking segments
%Produces DATAseg.mat files saved in subject-specific folders.
%By default, we have provided the walking segments that reproduce the 
%results of the paper.  This code will not overwrite existing DATAseg files.
%To rerun choose_walking_segments on your selected set, delete the 
%associated DATAseg files. This must be done manually.
choose_walking_segments(numericdata,metadata,0,redcap_path,proc_path,filesep);

%% Rotate data to earth-frame using sensor fusion
%Produces EARTHFRAME files saved in subject-specific folders. This code 
%will not overwrite existing EARTHFRAME files. To rerun sensor_earth_orient
%on your selected set, delete the associated EARTHFRAME files. This may be
%done manually or using clean_results.m
sensor_earth_orient(numericdata,metadata,redcap_path,proc_path,Vbl)

%% Integrate acceleration data
%Produces M2 files saved in subject-specific folders.  This code will not
%overwrite existing M2 files.  To rerun position_integration on your
%selected set, delete the associated M2 files.  This may be done manually
%or using clean_results.m
[pos, stridestats] = position_integration(numericdata,metadata,redcap_path,proc_path,Vbl);

% %% Compute and plot summary step length, duration, and speed
% compute_stride_metrics(metadata,redcap_path,proc_path,RCdata,RCind,record_id_vec,Vbl)

%% Reproduction of results and figures from paper
%In order to run this last section, the preceding code must be run to
%completion for all subjects for both lab and home recordings.  This
%reproduces the results and figures presented in the paper.
text = paper_results(redcap_path,proc_path,Vbl)