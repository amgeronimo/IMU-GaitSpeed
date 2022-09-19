%Uses the imufilter function of the MATLAB sensor fusion toolbox to rotate
%the data from sensor frame to earth frame.  Saves rotated data in
%EARTHFRAME files.
function [] = sensor_earth_orient(numericdata,metadata,redcap_path,proc_path,Vbl)

figureson = Vbl.figuresOn;
accHeader = Vbl.accHeader;
gyrHeader = Vbl.gyrHeader;
g_conversion = Vbl.g_conversion;

for i = 1:length(numericdata)
    pathPathfigFol = [redcap_path.bigOrg,filesep,redcap_path.indiv,filesep,metadata(i).SCname,filesep,metadata(i).SCname,proc_path.fig,filesep,metadata(i).SCname,proc_path.pathfig,filesep];
    pathDatafigFol = [redcap_path.bigOrg,filesep,redcap_path.indiv,filesep,metadata(i).SCname,filesep,metadata(i).SCname,proc_path.fig,filesep,metadata(i).SCname,proc_path.datafig,filesep];
    imuFol = [redcap_path.bigOrg,filesep,redcap_path.indiv,filesep,metadata(i).SCname,filesep,metadata(i).SCname,proc_path.imu,filesep];
    if sum(strcmp(metadata(i).datatype,{'Linear Acceleration','Accelerometer',accHeader}))>0
        disp(['File ' metadata(i).session]);

        %Find corresponding GYR data
        sidx = find(cellfun(@(x) strcmp(x,metadata(i).session),{metadata.session}));
        j = sidx(strcmp(gyrHeader,{metadata(sidx).datatype}));


        ai = []; qi = []; gi = [];
        tmp = cell(2,1);
        stationary = [];
        if exist([imuFol metadata(i).name '_EARTHFRAMEimu' '_' metadata(i).date '.mat'])>0
            disp([metadata(i).name '_EARTHFRAMEimu' '_' metadata(i).date '.mat <strong>exists</strong>']);
            load([imuFol metadata(i).name '_EARTHFRAMEimu' '_' metadata(i).date '.mat'])
            ai = accEO; qi = quat; gi = gyrEO;
        else
            clear walkseg
            dataseg = load([imuFol, metadata(i).name, '_DATAseg_',metadata(i).date,'.mat']);
            time = numericdata(i).time_interp_shift;
            acc = numericdata(i).pts_interp_shift/9.81;
            gyr = numericdata(j).pts_interp_shift*180/pi;
            slbl = {'acc','gyro'};
            dd_id = metadata(i).session;
            idx = tmp;
            stationary = stationary;
            walkseg = dataseg.walkseg;
            samplePeriod = time(2)-time(1);


            indexSel = logical(ones(length(time),1));
            time = time(indexSel);
            if ~isempty(gyr)
                gyrX = gyr(indexSel,1);
                gyrY = gyr(indexSel,2);
                gyrZ = gyr(indexSel,3);
            end
            if  ~isempty(acc)
                accX = acc(indexSel,1);
                accY = acc(indexSel,2);
                accZ = acc(indexSel,3);
            end


            % Acc and Gyro figure
            if figureson == 0
                signal_fig = figure('visible','off');
            else
                signal_fig = figure;
            end
            plot(time,acc); hold on;
            plot(time,gyr)
            xlabel('Time (s)')
            title(['Patient ' dd_id]);
            legend('ax','ay','az','gx','gy','gz')

            %% Detect stationary periods
            % Compute accelerometer magnitude
            acc_mag = sqrt(accX.*accX + accY.*accY + accZ.*accZ);

            % find value to subtract from acceleration data using the mode
            acc_center = mode(round(1000*acc_mag)/1000);
            % acc_stationary = mean(acc_center);

            % subtracting mean acceleration from acceleration magnitude
            acc_magD = acc_mag-mean(acc_center); % check why mean (mode)

            % HP filter accelerometer data
            filtCutOff = 0.0001;
            [b, a] = butter(1, (2*filtCutOff)/(1/samplePeriod), 'high');
            acc_magFilt_h = filtfilt(b, a, [zeros(1000,1); acc_magD; zeros(1000,1)]);

            % Compute absolute value
            acc_magFilt_h = abs(acc_magFilt_h(1001:end-1000));

            % LP filter accelerometer data
            filtCutOff = 5;
            [b, a] = butter(1, (2*filtCutOff)/(1/samplePeriod), 'low');
            acc_magFilt = filtfilt(b, a, acc_magFilt_h);

            % find peaks filter on acc x, y, z data
            filtCutOff_FP = 2;
            [b, a] = butter(1, (2*filtCutOff)/(1/samplePeriod), 'low');
            accX_Filt = filtfilt(b, a, accX);
            accY_Filt = filtfilt(b, a, accY);
            accZ_Filt = filtfilt(b, a, accZ);

            % acc Filt figure
            if figureson == 0
                accFilt_fig = figure('visible','off');
            else
                accFilt_fig = figure;
            end
            plot(time,acc_magFilt);
            hold on
            plot(time,accX_Filt);
            plot(time,accY_Filt);
            plot(time,accZ_Filt);
            hold off
            title(['Patient ' dd_id])
            xlabel('Time (s)')
            ylabel('Acc Filt')
            legend('acc mag Filt','acc X Filt','acc Y Filt','acc Z Filt');


            %Remove bad walking segments (less than 5s)
            % start ws_mult at 4 and change if the criteria is met below
            ws_mult = 4;
            if dd_id(1) == 'C' && dd_id(4) == 't' % if control and gait and motor file
                % do nothing to walkseg
                % set multiplier for sd comparison to 2
                ws_mult = 2;
            else
                if dd_id(1) == 'S' && dd_id(5) == 't' % if patient and gait and motor file
                    % find bad walking segments less than 5 seconds and remove
                    badind = diff(walkseg,[],2)<=250;
                    walkseg = walkseg(~badind,:);
                    % set multiplier for sd comparison to 4
                    ws_mult = 4;
                end
            end

            % find peaks on (+) and (-) acc xyz data AND acc_magFilt
            locs_posX = cell(size(walkseg,1),1);
            locs_negX = cell(size(walkseg,1),1);
            locs_posY = cell(size(walkseg,1),1);
            locs_negY = cell(size(walkseg,1),1);
            locs_posZ = cell(size(walkseg,1),1);
            locs_negZ = cell(size(walkseg,1),1);
            locs_mag = cell(size(walkseg,1),1);

            % determine if the subject is sc8 to set the min peak distance to 50
            if dd_id(regexpi(dd_id(1:4),'[^SC-]')) == '8'
                mpd = 50;
            else
                mpd = 35;
            end

            % extract the acc for only the walking segments and find peaks
            for k = 1:size(walkseg,1)
                ws = sort(walkseg(k,:));
                accX_pos = accX_Filt(ws(1):ws(2));
                accX_neg = -accX_Filt(ws(1):ws(2));
                accY_pos = accY_Filt(ws(1):ws(2));
                accY_neg = -accY_Filt(ws(1):ws(2));
                accZ_pos = accZ_Filt(ws(1):ws(2));
                accZ_neg = -accZ_Filt(ws(1):ws(2));
                acc_magFilt_FP = acc_magFilt(ws(1):ws(2));
                [~,locs_posX{k}] = findpeaks(accX_pos,'MinPeakHeight',mean(accX_pos)+std(accX_pos),'MinPeakDistance',mpd);
                [~,locs_negX{k}] = findpeaks(accX_neg,'MinPeakHeight',mean(accX_neg)+std(accX_neg),'MinPeakDistance',mpd);
                [~,locs_posY{k}] = findpeaks(accY_pos,'MinPeakHeight',mean(accY_pos)+std(accY_pos),'MinPeakDistance',mpd);
                [~,locs_negY{k}] = findpeaks(accY_neg,'MinPeakHeight',mean(accY_neg)+std(accY_neg),'MinPeakDistance',mpd);
                [~,locs_posZ{k}] = findpeaks(accZ_pos,'MinPeakHeight',mean(accZ_pos)+std(accZ_pos),'MinPeakDistance',mpd);
                [~,locs_negZ{k}] = findpeaks(accZ_neg,'MinPeakHeight',mean(accZ_neg)+std(accZ_neg),'MinPeakDistance',mpd);
                [~,locs_mag{k}] = findpeaks(acc_magFilt_FP,'MinPeakHeight',mean(acc_magFilt_FP)+std(acc_magFilt_FP),'MinPeakDistance',mpd);
            end

            % calculate step duration based on findpeaks and determine sd
            diff_posX = cellfun(@diff,locs_posX,'UniformOutput',false);
            diff_negX = cellfun(@diff,locs_negX,'UniformOutput',false);
            diff_posY = cellfun(@diff,locs_posY,'UniformOutput',false);
            diff_negY = cellfun(@diff,locs_negY,'UniformOutput',false);
            diff_posZ = cellfun(@diff,locs_posZ,'UniformOutput',false);
            diff_negZ = cellfun(@diff,locs_negZ,'UniformOutput',false);
            diff_mag = cellfun(@diff,locs_mag,'UniformOutput',false);
            dur_posX = cell2mat(diff_posX)*0.02;
            dur_negX = cell2mat(diff_negX)*0.02;
            dur_posY = cell2mat(diff_posY)*0.02;
            dur_negY = cell2mat(diff_negY)*0.02;
            dur_posZ = cell2mat(diff_posZ)*0.02;
            dur_negZ = cell2mat(diff_negZ)*0.02;
            dur_mag = cell2mat(diff_mag)*0.02;

            % this will make all sd's equal to Nan's for gait and motor assessment
            % files especially for controls. need to adjust the multiplier to 2 or 3 for
            % controls, and leave as 4 for patients
            % there should be at least 2 steps per walking segment. there's always 3 ws
            % for gm files.
            if length(dur_posX)<(ws_mult*size(walkseg,1)) % why is this one 2
                sd_posX = NaN;
            else
                sd_posX = std(dur_posX);
            end

            if length(dur_negX)<(ws_mult*size(walkseg,1))
                sd_negX = NaN;
            else
                sd_negX = std(dur_negX);
            end
            if length(dur_posY)<(ws_mult*size(walkseg,1))
                sd_posY = NaN;
            else
                sd_posY = std(dur_posY);
            end
            if length(dur_negY)<(ws_mult*size(walkseg,1))
                sd_negY = NaN;
            else
                sd_negY = std(dur_negY);
            end
            if length(dur_posZ)<(ws_mult*size(walkseg,1))
                sd_posZ = NaN;
            else
                sd_posZ = std(dur_posZ);
            end
            if length(dur_negZ)<(ws_mult*size(walkseg,1))
                sd_negZ = NaN;
            else
                sd_negZ = std(dur_negZ);
            end
            if length(dur_mag)<(ws_mult*size(walkseg,1))
                sd_mag = NaN;
            else
                sd_mag = std(dur_mag);
            end

            % determine which sd is the lowest
            sd_all = [sd_posX;sd_negX;sd_posY;sd_negY;sd_posZ;sd_negZ;sd_mag];
            sd_all(sd_all==0) = NaN; % set sd==0 to NaN
            [~,minind] = min(sd_all,[],'omitnan'); % ignore NaNs
            % determine what set of acc and peaks to use
            if minind == 1 % posx had the lowest sd
                acc_Filt_FP = accX_Filt;
                locs = locs_posX;
                acc_id = 'accXpos';
            else
                if minind == 2 % negx had the lowest sd
                    acc_Filt_FP = -accX_Filt;
                    locs = locs_negX;
                    acc_id = 'accXneg';
                else
                    if minind == 3 % posy had the lowest sd
                        acc_Filt_FP = accY_Filt;
                        locs = locs_posY;
                        acc_id = 'accYpos';
                    else
                        if minind == 4 % negy had the lowest sd
                            acc_Filt_FP = -accY_Filt;
                            locs = locs_negY;
                            acc_id = 'accYneg';
                        else
                            if minind == 5 % posz had the lowest sd
                                acc_Filt_FP = accZ_Filt;
                                locs = locs_posZ;
                                acc_id = 'accZpos';
                            else
                                if minind == 6 % negz had the lowest sd
                                    acc_Filt_FP = -accZ_Filt;
                                    locs = locs_negZ;
                                    acc_id = 'accZneg';
                                else
                                    if minind == 7 % accmag had the lowest sd
                                        acc_Filt_FP = acc_magFilt;
                                        locs = locs_mag;
                                        acc_id = 'accmag';
                                    end
                                end
                            end
                        end
                    end
                end
            end

            % plot the peaks on the walksegs and trim the locs
            if figureson == 0
                peaks_fig = figure('visible','off');
            else
                peaks_fig = figure;
            end
            locs_trim = cell(size(walkseg,1),1);
            plot(time,acc_Filt_FP)
            hold on
            for k = 1:size(walkseg,1)
                plot(time(walkseg(k)+locs{k}-1),acc_Filt_FP(walkseg(k)+locs{k}-1),'r*')
                locs_trim{k} = locs{k}(2:end-1);
            end
            hold off
            title(['Find Peaks on filtered ' acc_id 'data ' dd_id])
            xlabel('Time (s)')
            ylabel(acc_id)


            % calculate peak to peak duration on trimmed locs
            diff_peaks = cellfun(@diff,locs_trim,'UniformOutput',false);
            diff_peaks(cellfun(@(x) sum(isempty(x))==1,diff_peaks)) = [];
            dur_peaks = cell2mat(diff_peaks)*0.02;
            % remove outliers here
            outliers_low = dur_peaks(dur_peaks>0.2);
            outliers_max = mean(outliers_low)+std(outliers_low);
            good_steps = outliers_low(outliers_low<outliers_max);
            dur_peaks_mean = mean(good_steps); % new method, removing outliers

            % plot histogram of peak durations
            if figureson == 0
                durhist_fig = figure('visible','off');
            else
                durhist_fig = figure;
            end
            histogram(dur_peaks,0:0.1:3)
            title(['Peak to Peak Stride Duration - ' dd_id])
            xlabel('Peak to Peak Duration (s)')

            % find peaks on entire signal
            [~,locs_all] = findpeaks(acc_Filt_FP,'MinPeakHeight',mean(acc_Filt_FP)+std(acc_Filt_FP),'MinPeakDistance',mpd);

            % find minimum in between peaks for entire signal
            if figureson == 0
                min_fig = figure('visible','off');
            else
                min_fig = figure;
            end
            plot(time,acc_magFilt)
            hold on
            stationary = zeros(length(acc_magFilt),1);
            minval = zeros(length(locs_all),1);
            minloc = zeros(length(locs_all),1);
            % loop through each loc (peak)
            for j = 1:length(locs_all)-1
                % find min in between peaks
                dist = locs_all(j):(locs_all(j)+(0.75*(locs_all(j+1)-locs_all(j))));
                [minval(j),minloc(j)] = min(acc_magFilt(dist));
                stationary(locs_all(j)+minloc(j)) = 1;
                plot(time(locs_all(j)+minloc(j)),minval(j),'b*')
            end
            hold off
            title(['Minimum in between peaks ' dd_id])
            xlabel('Time (s)')
            ylabel(acc_id)


            % plot threshold onto acc
            if figureson == 0
                thresh_fig = figure('visible','off');
            else
                thresh_fig = figure;
            end
            plot(time, acc_magFilt, 'r');
            hold on
            ss=plot(time, .2*stationary, 'k');
            title('Acc data with selected threshold')
            xlabel('Time (s)');
            ylabel('Acceleration (g)');
            ss.YData = .2*double(stationary);
            ws_time = round(0.02*walkseg);
            for j = 1:size(ws_time,1)
                xws = [ws_time(j,1) ws_time(j,2) ws_time(j,2) ws_time(j,1)];
                yws = [0 0 max(acc_magFilt) max(acc_magFilt)];
                patch(xws,yws,'k','FaceAlpha',0.1,'EdgeColor','none')
            end
            hold off
            title(['Threshold ' dd_id])
            xlabel('Time (s)')
            ylabel('acc mag Filt')

            % peaks and min plot
            peaks_mins = figure;
            subplot(2,1,1)
            plot(time,acc_Filt_FP)
            hold on
            plot(time(locs_all),acc_Filt_FP(locs_all),'r*')
            hold off
            title(acc_id)
            subplot(2,1,2)
            plot(time, acc_magFilt)
            hold on
            plot(time(locs_all+minloc),minval,'b*')
            hold off
            title('acc mag filt')
            linkaxes

            lin_acc_noise =  0.0096236;%(m/s²)²
            lin_acc_decay = 0.5;

            %Values from data sheets
            acc_noise = 1.8*g_conversion/10^3; %From Bosch BMI160 datasheet 1.8 mg rms;
            gyro_noise = (deg2rad(0.07))^2; %From Bosch BMI160 datasheet


            clear ifilt
            ifilt = imufilter('SampleRate', 1/samplePeriod,'AccelerometerNoise',acc_noise,...
                'GyroscopeNoise',gyro_noise,...
                'LinearAccelerationNoise',lin_acc_noise*100,'LinearAccelerationDecayFactor',lin_acc_decay);
            [qimu, av] = ifilt([accX accY accZ]*g_conversion, deg2rad([gyrX gyrY gyrZ]));
            quat = compact(qimu);

            % Compute translational accelerations
            %matlab built-in takes quat --- this rotation is into NED coordinates
            accEO=rotatepoint(quaternion(quat),[accX accY accZ]);
            gyrEO = rad2deg(av);
            timeEO = time;

            % Convert acceleration measurements to m/s/s
            accEO = accEO * g_conversion;

            hh=figure('Position',[0 0 600 600]);
            plot_position(dd_id,accEO,samplePeriod,stationary,walkseg);
            sgtitle(dd_id);

            % save figures
            saveas(signal_fig,[pathDatafigFol metadata(i).name '_AccGyr_' metadata(i).date '.pdf'])
            saveas(accFilt_fig,[pathDatafigFol metadata(i).name '_accFilt_' metadata(i).date '.fig'])
            saveas(peaks_fig,[pathDatafigFol metadata(i).name '_peaks_' metadata(i).date '.fig'])
            saveas(durhist_fig,[pathDatafigFol metadata(i).name '_durHist_' metadata(i).date '.pdf'])
            saveas(min_fig,[pathDatafigFol metadata(i).name '_minbetweenpeaks_' metadata(i).date '.fig'])
            saveas(thresh_fig,[pathDatafigFol metadata(i).name '_threshold_' metadata(i).date '.fig'])
            saveas(peaks_mins,[pathDatafigFol metadata(i).name 'peaks_mins_' metadata(i).date '.fig'])
            saveas(hh,[pathPathfigFol metadata(i).name '_ORTHOimu_' metadata(i).date '.fig'])
            metadataEO = metadata(i);

            save([imuFol metadata(i).name '_EARTHFRAMEimu' '_' metadata(i).date '.mat'],...
                'accEO','gyrEO','timeEO','quat','metadataEO','stationary','acc_id');

        end
        if figureson == 1
            gg = figure('Visible','off','Position',[0 100 1200 400]);
            plot(accEO); hold on;
            plot(numericdata(i).pts_interp_shift);
            legend('AP','ML','VV','X','Y','Z')
            load([imuFol, metadataEO.name '_DATAseg_' metadata(i).date])
            if ~isempty(walkseg)
                if walkseg(1,1) < walkseg(1,2)
                    xlim([walkseg(1,1) walkseg(1,2)]);
                else
                    xlim([walkseg(1,2) walkseg(1,1)]);
                end
            end
            saveas(gcf,[pathPathfigFol metadata(i).name '_EARTHFRAMEimu_' metadata(i).date '.pdf']);
        end
        close all;
    end
end