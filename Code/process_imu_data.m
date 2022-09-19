%Data is saved on per-subject bases in SCXX_ALLDATA.mat
% numericdata contains timestamp vectors for recordings
% metadata contains files' names, dates, assistance code, etc.

function [data, metadata] = process_imu_data(redcap_path,proc_path,Vbl)
disp '     Begin processing IMU data';

subj = Vbl.subj;
sctmp = regexpi(subj,'[a-z]+','match');
sntmp = regexp(subj,'\d+','match');
if length(sctmp{1})==2 %Contains TMW suffix
    savesubj = cellfun(@(x,y,z) [upper(x) sprintf('%02d',str2double(y)) z],...
        cellfun(@(x) x(1),sctmp),[sntmp{:}],  cellfun(@(x) x(2),sctmp),'UniformOutput',false);
else
    savesubj = cellfun(@(x,y) [upper(x) sprintf('%02d',str2double(y))],[sctmp{:}],[sntmp{:}],'UniformOutput',false);
end
figureson = Vbl.figuresOn;
path = Vbl.pathtodata;
acc = Vbl.accHeader;
gyr = Vbl.gyrHeader;

warning('off');

bigOrg = redcap_path.bigOrg;
rawall = redcap_path.rawall;
indiv = redcap_path.indiv;
path_fig = proc_path.fig;
path_preprocfig = proc_path.preprocfig;
path_csv = proc_path.csv;

%% Load Metadata
for si = 1:length(savesubj)
    clear metadata data OD
    j = 1;
    k=1;
    tmpcomp = [];
    Mfiles = dir([path subj{si} '_*']);
    Dfile = dir([path savesubj{si} 'DATA.mat']);
 
    if isempty(Dfile)
        disp([[savesubj{si} 'DATA.mat'] ' is not up to date or does not exist']);
        OD.names.metadata.filename = '';
    else
        OD.names = load([path savesubj{si} 'DATA.mat']);
    end

    for i = 1:length(Mfiles)
        tmpfile = Mfiles(i);

        if strcmp(tmpfile.name(end),'v')                                        % csv files, ignoring '.' and '..'
            tmpi = strfind(tmpfile.name,'_');
            if ~isempty(tmpi)
                metadatatmp.name = tmpfile.name(1:tmpi(1)-1);                   % patient's whole study code
                metadatatmp.date = tmpfile.name(tmpi(2)+1:tmpi(3)-1);           % contains date from new file naming convention
                metadatatmp.datatype = tmpfile.name(tmpi(3)+1:tmpi(3)+3);       % GYR or ACC
                if dir([path,filesep,tmpfile.name]).bytes < 10000
                    disp(['Ignoring ',tmpfile.name,', consider deleting session']);
                elseif sum(strcmp(tmpcomp,[metadatatmp.name,metadatatmp.date(1:8),metadatatmp.datatype,num2str(tmpfile.bytes)]))>0
                    % if after first file, check if study code, date, file type, and exact file size
                    % are duplicated. --HL 1/27/2021
                    duplicatefiles{k} = tmpfile.name;
                    disp(['Ignoring duplicate file ',tmpfile.name,', consider deleting session']);
                    k=k+1;
                else
                    metadata(j).filename = tmpfile.name;                            % whole file name
                    metadata(j).name = tmpfile.name(1:tmpi(1)-1);                   % patient's whole study code
                    tmpt = regexp(metadata(j).name,'[a-zA-Z]*','match');
                    tmpn = regexp(metadata(j).name,'\d*','match');
                    if length(tmpt)>1 %files ending in tmw
                        metadata(j).SCname = metadata(j).name;
                    else
                        metadata(j).SCname = sprintf('%s%02d',upper(tmpt{:}),str2num(cell2mat(tmpn)));
                    end

                    metadata(j).date = tmpfile.name(tmpi(2)+1:tmpi(3)-1);           % contains date from new file naming convention
                    metadata(j).assist = tmpfile.name(tmpi(1)+1:tmpi(2)-1);         % contains mode of assistance from new file naming convention
                    metadata(j).session = [metadata(j).name,'-',metadata(j).date];  % unique series of #'s corresponding to this file's session
                    metadata(j).datatype = tmpfile.name(tmpi(3)+1:tmpi(3)+3);       % GYR or ACC
                    % add field for iPhone (15)/Android (18)
                    if length(metadata(j).date) > 15 % android
                        metadata(j).phonetype = 'Android';
                    else %iPhone
                        metadata(j).phonetype = 'iPhone';
                    end
                    tmpcomp{j} = [metadata(j).name,metadata(j).date(1:8),metadata(j).datatype,num2str(tmpfile.bytes)];
                    j=j+1;
                end
            end
        end
    end
    clear tmpdate tmpfile;

    if ~isempty(Dfile)
        if length(OD.names.metadata)==length(metadata)
            disp([[savesubj{si} 'DATA.mat'] ' up to date']);
            data = OD.names.data;
            metadata = OD.names.metadata;
            continue;
        end
    end


    %% Load Data
    %Interpolate to a uniform sampling rate (50 Hz) and shift data so that
    %the times of each sensor of the two devices match up, trim to equal
    %length.
    j = 1;
    for i = j:length(metadata)
        idx = contains({OD.names.metadata.filename},metadata(i).filename);
        if sum(idx) == 0
            if isempty(OD.names.metadata.filename)
                disp([num2str(j) '/' num2str(length(metadata)) '  (' metadata(i).filename ')'])
            else
                disp([num2str(j) '/' num2str(length(metadata) - size({OD.names.metadata.filename},2))])
            end
            [~,~,raw1] = xlsread([path metadata(i).filename]);
            dataind = find(contains(raw1(1,:),'x','IgnoreCase',true),1);
            datavec = cell2mat(raw1(2:end,dataind:dataind+2));
            switch metadata(i).datatype
                case acc
                    data(i).pts = datavec*Vbl.g_conversion;
                case gyr
                    data(i).pts = datavec*pi/180;
            end
            timevec = raw1(2:end,1);

            %Date conversions-- iOS vs Android
            if contains(num2str(cell2mat(timevec(1))),'-')==true
                tmptimenum = cellfun(@(x) [x(1:10) '.' x(12:23)],timevec,'UniformOutput',false);
                data(i).timenum = cellfun(@(x) datenum(x,'yyyy-mm-dd.HH:MM:SS.FFF'),tmptimenum);
            else
                tmp = datetime(cell2mat(timevec)/1000,'ConvertFrom','posixtime','Timezone','EST','format','yyyyMMdd-HH:mm:ss:SSS');
                data(i).timenum = datenum(tmp);
            end

            data(i).starttime = data(i).timenum(1);
            data(i).elaptime = (data(i).timenum-data(i).timenum(1))*3600*24;

            obsrate = 1/mean(diff(data(i).elaptime));
            expectedrate = 50;
            if abs(obsrate - expectedrate) >= 1
                disp(['The observed rate is ' num2str(obsrate),...
                    'Hz for the ' metadata(i).datatype ' file of recording ' metadata(i).session]);
            end
            metadata(i).fs = expectedrate;

            %% Set even spacing between samples
            % Need to adjust the time within samples to be evenly spaced

            %Resample all data to defined rate
            fsi = 50;
            metadata(i).fs_interp = fsi;

            %Resample data to even time spacing
            data(i).time_uniform = (1:length(data(i).pts))/obsrate;
            [data(i).pts_interp, data(i).time_interp] = resample(data(i).pts,...
                data(i).time_uniform,metadata(i).fs_interp);

            %% Enforce same start time and recording length across sensors
            tmp = find(strcmp({metadata.session},metadata(i).session));
            if i==tmp(end)

                %Find latest starttime and number of samples to remove from
                %beginning of other files
                sts = cat(1,data(tmp).starttime);
                sts_rem =  round(str2num(datestr(max(sts)-sts,'SS.FFF'))*metadata(i).fs_interp);
                kk = 1;
                for ix = 1:length(tmp)
                    data(tmp(ix)).pts_interp_shift = data(tmp(ix)).pts_interp(sts_rem(kk)+1:end,:);
                    data(tmp(ix)).time_interp_shift = data(tmp(ix)).time_interp(sts_rem(kk)+1:end)-data(tmp(ix)).time_interp(sts_rem(kk)+1);
                    data(tmp(ix)).starttime = max(sts);
                    kk=kk+1;
                end

                %Find shortest file and trim all to that length
                tmpend = min(cellfun(@length,{data(tmp).time_interp_shift}));
                for ix = 1:length(tmp)
                    data(tmp(ix)).pts_interp_shift = data(tmp(ix)).pts_interp_shift(1:tmpend,:);
                    data(tmp(ix)).time_interp_shift = data(tmp(ix)).time_interp_shift(1:tmpend);
                    % copies pts_interp for acc and gyr to calibrated
                    data(tmp(ix)).pts_calibrated = data(tmp(ix)).pts_interp_shift;
                    %Calibrate Magnetometer data


                end

                if figureson == 0
                    imu_fig = figure('visible','off');
                else
                    imu_fig = figure;
                end

                xxi = 1;
                for ix = 1:length(tmp)-1
                    down = 3;
                    xxxi(xxi)=subplot(length(tmp)-1,1,xxi); hold on;
                    title(['Patient ' metadata(tmp(ix)).name ' ' metadata(tmp(ix)).datatype]);

                    plot(data(tmp(ix)).elaptime,...
                        data(tmp(ix)).pts(:,down),'r');
                    plot(data(tmp(ix)).time_interp,...
                        data(tmp(ix)).pts_interp(:,down),'b');
                    plot(data(tmp(ix)).time_interp_shift,...
                        data(tmp(ix)).pts_interp_shift(:,down),'Color',[0 .8 0]);
                    xxi=xxi+1;
                end
                legend('raw data',[num2str(fsi) ' Hz resample'],[num2str(fsi) ' Hz resample w/time shift']);
                linkaxes(xxxi,'x');

                tmp = get(imu_fig,'Position');
                set(imu_fig,...
                    'DefaultAxesFontSize',5,...
                    'DefaultLineLineWidth',2,...
                    'PaperUnits','points',...
                    'PaperSize',[tmp(3) tmp(4)],... 
                    'PaperPosition',[0 0 tmp(3) tmp(4)]); 
                saveas(imu_fig,[bigOrg,filesep,indiv,filesep,metadata(i).SCname,filesep,metadata(i).SCname,path_fig,filesep,metadata(i).SCname,path_preprocfig,filesep,metadata(i).SCname '_proc_' metadata(i).date '.pdf'],'pdf');
                close(imu_fig);
            else
                data(i).time_interp_shift = data(i).time_interp;
                data(i).pts_interp_shift = data(i).pts_interp;
            end

        else
            ODidx = find(idx);
            data(i) = OD.names.data(ODidx);
        end
        j = j + 1;
    end
    save([path [savesubj{si} 'DATA.mat']],'metadata','data')
end

%% Pass all data out of the function (may be data from multiple subjects)
data = []; metadata = [];
for si = 1:length(savesubj)
    SD = load([path savesubj{si} 'DATA.mat']);
    data = cat(2,data,SD.data);
    metadata = cat(2,metadata,SD.metadata);
end