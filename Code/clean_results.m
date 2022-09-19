%Clean folder structure

rec_i = input('Clean results of home (0) or lab (1) recordings?\n');

if rec_i==0
    pathtodata = ['..' filesep 'HomeRecordings' filesep];

elseif rec_i==1
    pathtodata = ['..' filesep 'LabClinicRecordings' filesep];
else
    disp('Invalid selection');
    return;
end
pathtoresults = 'Processed_Data\Individual_Patient_Data';


x = input('Clean the results of process_imu_data ("DATA" files)? yes (1) or no (0)\n');
if x

    F = dir([pathtodata '*DATA.mat']);
    Fnames = {F.name};
    Snames = regexpi(Fnames,'s?c\d+.*(?=DATA.mat)','match');
    for i = 1:length(Snames)
        fprintf(['(' num2str(i) ') ' Snames{i}{:} '\n'])
    end
    subj_i = input('Type the number next to the subject of the DATA file to delete.  Type 0 for all. \n');
    if subj_i==0
        subj = [Snames{:}];
    elseif subj_i<=length(Snames)
        subj = Snames{subj_i};
    else
        disp('Invalid selection');
        return;
    end

    if x
        for si = 1:length(subj)
            if exist([pathtodata subj{si} 'DATA.mat'])==2
                delete([pathtodata subj{si} 'DATA.mat'])
            end
        end
    end
end

%Remove DATAseg walking segments?
%Not completing because we want these to be consistent

%Remove EARTHORIENT files?
x = input('Clean results of sensor_earth_orient ("EARTHFRAME" files)? yes (1) or no (0)?\n');
if x


    F = dir(pathtoresults);
    Fnames = {F.name};
    Snames = regexpi(Fnames,'s?c\d+(tmw[LR])*','match');
    Snames = [Snames{:}]
    for i = 1:length(Snames)
        fprintf(['(' num2str(i) ') ' Snames{i} '\n'])
    end
    subj_i = input('Type the number next to the subject of the EARTHFRAME files to delete.  Type 0 for all. \n');
    if subj_i==0
        subj = Snames;
    elseif subj_i<=length(Snames)
        subj = Snames(subj_i);
    else
        disp('Invalid selection');
        return;
    end



    for si = 1:length(subj)
        ftmp = dir([pathtoresults filesep subj{si} filesep '**\*EARTHFRAME*']);
        for fi = 1:length(ftmp)
            delete([ftmp(fi).folder filesep ftmp(fi).name])
        end
    end
end

%Remove M2 files?
x = input('Clean results of position_integration ("M2" files)? yes (1) or no (0)?\n');
if x

    F = dir(pathtoresults);
    Fnames = {F.name};
    Snames = regexpi(Fnames,'s?c\d+(tmw[LR])*','match');
    Snames = [Snames{:}]
    for i = 1:length(Snames)
        fprintf(['(' num2str(i) ') ' Snames{i} '\n'])
    end
    subj_i = input('Type the number next to the subject of the M2 files to delete.  Type 0 for all. \n');
    if subj_i==0
        subj = Snames;
    elseif subj_i<=length(Snames)
        subj = Snames(subj_i);
    else
        disp('Invalid selection');
        return;
    end

    for si = 1:length(subj)
        ftmp = dir([pathtoresults filesep subj{si} filesep '**\' 'M2parameters*']);
        for fi = 1:length(ftmp)
            delete([ftmp(fi).folder filesep ftmp(fi).name])
        end
    end
end

%Remove Figures?
x = input('Clean figure files? yes (1) or no (0)?\n');
if x

    F = dir(pathtoresults);
    Fnames = {F.name};
    Snames = regexpi(Fnames,'s?c\d+(tmw[LR])*','match');
    Snames = [Snames{:}]
    for i = 1:length(Snames)
        fprintf(['(' num2str(i) ') ' Snames{i} '\n'])
    end
    subj_i = input('Type the number next to the subject of the figure files to delete.  Type 0 for all. \n');
    if subj_i==0
        subj = Snames;
    elseif subj_i<=length(Snames)
        subj = Snames(subj_i);
    else
        disp('Invalid selection');
        return;
    end

    for si = 1:length(subj)
        ftmp = dir([pathtoresults filesep subj{si} filesep subj{si} '_Figures\**\*.fig']);
        for fi = 1:length(ftmp)
            delete([ftmp(fi).folder filesep ftmp(fi).name])
        end
        ftmp = dir([pathtoresults filesep subj{si} filesep subj{si} '_Figures\**\*.pdf']);
        for fi = 1:length(ftmp)
            delete([ftmp(fi).folder filesep ftmp(fi).name])
        end
    end
end
