%This function integrates the rotated acceration to velocity.  Zero
%velocity update is performed at each point the foot is determined to be
%stationary.  Position is integrated from velocity, from which individual
%stride lengths, durations, and speeds are calculated.
function [posPlot,stridestats] = position_integration(data,metadata,redcap_path,proc_path,Vbl)

figureson = Vbl.figuresOn;
accHeader = Vbl.accHeader;
gyrHeader = Vbl.gyrHeader;

bigOrg = redcap_path.bigOrg;
rawall = redcap_path.rawall;
indiv = redcap_path.indiv;
group = redcap_path.group;

path_csv = proc_path.csv;
path_fig = proc_path.fig;
path_results = proc_path.results;
path_datafig = proc_path.datafig;
path_pathfig = proc_path.pathfig;
path_imu = proc_path.imu;
path_step = proc_path.step;
path_line = proc_path.line;

g_conversion = Vbl.g_conversion;
earthframe = 1;

for i = 1:length(data)
    disp([num2str(i) '/' num2str(length(data))])
    close all;
    if strcmp(metadata(i).datatype,accHeader)
        pathfigFol = [redcap_path.bigOrg,filesep,redcap_path.indiv,filesep,metadata(i).SCname,filesep,metadata(i).SCname,path_fig,filesep,metadata(i).SCname,path_pathfig,filesep];
        imuFol = [redcap_path.bigOrg,filesep,redcap_path.indiv,filesep,metadata(i).SCname,filesep,metadata(i).SCname,path_imu,filesep];
        stepFol = [redcap_path.bigOrg,filesep,redcap_path.indiv,filesep,metadata(i).SCname,filesep,metadata(i).SCname,path_step,filesep];
        if earthframe==1
            try
                load([imuFol metadata(i).name '_EARTHFRAMEimu' '_' metadata(i).date '.mat'])
            catch
                disp(['No EARTHFRAME data for ' metadata(i).session]);
                continue;
            end
            acc = accEO;
            time = timeEO;
        else
            error('Needs earth frame data');
        end
        samplePeriod = time(2)-time(1);

        saveL = ['imu_' metadata(i).session];

        if exist ([stepFol 'M2parameters_' saveL '.mat'])>0
            M2parameters = load([stepFol 'M2parameters_' saveL '.mat']);
            posPlot = M2parameters.posPlot;
            stationary = M2parameters.stationary;
            stridestats = M2parameters.stridestats;
            disp (['M2parameters file <strong>exists</strong> for ' metadata(i).session]);
        else
            % -------------------------------------------------------------------------
            % Integrate acceleration to yield velocity
            acc(:,3) = acc(:,3) - g_conversion;
            vel = zeros(size(acc));
            for t = 2:length(vel)
                vel(t,:) = vel(t-1,:) + acc(t,:) * samplePeriod;
                if(stationary(t) == 1)
                    vel(t,:) = [0 0 0];     % force zero velocity when foot stationary
                end
            end

            % Compute integral drift during non-stationary periods
            velDrift = zeros(size(vel));
            stationaryStart = find([0; diff(stationary)] == -1);
            stationaryEnd = find([0; diff(stationary)] == 1);
            if length(stationaryEnd)>length(stationaryStart) % ag 7/30/19
                stationaryEnd(1) = [];
            end
            stationaryStart = stationaryStart(1:length(stationaryEnd)); %ag 5/9/19
            if stationaryEnd(1)<stationaryStart(1)
                stationaryEnd = stationaryEnd(2:end);
                stationaryStart = stationaryStart(1:end-1);
                for si = 1:numel(stationaryEnd)
                    driftRate = vel(stationaryEnd(si)-1, :) / (stationaryEnd(si) - stationaryStart(si));
                    enum = 1:(stationaryEnd(si) - stationaryStart(si));
                    drift = [enum'*driftRate(1) enum'*driftRate(2) enum'*driftRate(3)];
                    velDrift(stationaryStart(si):stationaryEnd(si)-1, :) = drift;
                end
            else
                for si = 1:numel(stationaryEnd)
                    driftRate = vel(stationaryEnd(si)-1, :) / (stationaryEnd(si) - stationaryStart(si));
                    enum = 1:(stationaryEnd(si) - stationaryStart(si));
                    drift = [enum'*driftRate(1) enum'*driftRate(2) enum'*driftRate(3)];
                    velDrift(stationaryStart(si):stationaryEnd(si)-1, :) = drift;
                end
            end

            % Remove integral drift
            vel = vel - velDrift;


            % -------------------------------------------------------------------------
            % Integrate velocity to yield position
            pos = zeros(size(vel));
            for t = 2:length(pos)
                pos(t,:) = pos(t-1,:) + vel(t,:) * samplePeriod;    % integrate velocity to yield position
            end

            % Plot translational position
            if figureson == 0
                pos_fig = figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', 'Position','visible','off');
            else
                pos_fig = figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', 'Position');
            end
            plot(time, pos(:,1), 'r');
            hold on;
            plot(time, pos(:,2), 'g');
            plot(time, pos(:,3), 'b');
            title(['Position: ',metadata(i).session]);
            xlabel('Time (s)');
            ylabel('Position (m)');
            legend('X', 'Y', 'Z');
            hold off;
            saveas(pos_fig,[pathfigFol metadata(i).name '_trans_position_' metadata(i).date '.pdf'])

            posPlot = pos;

            %Load data periods of walking
            clear walkseg
            dataseg = load([imuFol, metadata(i).name, '_DATAseg_',metadata(i).date,'.mat']);
            walkseg = dataseg.walkseg;

            %Remove bad walking segments (less than 5s) from home recordings
            findtmw = strfind(stepFol,'tmw');
            if isempty(findtmw) % if NOT a gait and motor file
                % find bad walking segments less than 5 seconds and remove
                badind = diff(walkseg,[],2)<=250;
                walkseg = walkseg(~badind,:);
                % otherwise, leave walkseg alone
            end




            kpp = zeros(length(stationaryStart),1);
            for wi = 1:size(walkseg,1)
                patch([walkseg(wi,1) walkseg(wi,2) walkseg(wi,2) walkseg(wi,1) walkseg(wi,1)],...
                    [-30 -30 30 30 -30],'k','EdgeAlpha',0,'FaceAlpha',.1)
                kpp = kpp | (stationaryStart>walkseg(wi,1) & stationaryEnd<walkseg(wi,2));
            end


            % overview path plot
            cmpb = zeros(size(posPlot));
            if figureson == 0
                ovpath_fig = figure('Position',[0 0 600 600],'renderer','painters','visible','off');
            else
                ovpath_fig = figure('Position',[0 0 600 600],'renderer','painters');
            end
            ox = scatter3(posPlot(:,1),posPlot(:,2),posPlot(:,3),10,cmpb,'Marker','.');
            hold on
            scatter3(posPlot(1,1),posPlot(1,2),posPlot(1,3),50,[0 0 0],'Marker','o');
            scatter3(posPlot(end,1),posPlot(end,2),posPlot(end,3),50,[0 0 0],'Marker','x');
            for wi = 1:size(walkseg,1)
                ox.CData(walkseg(wi,1):walkseg(wi,2),:) = repmat([1 0 0],diff(walkseg(wi,:))+1,1);
            end
            xl = get(gca,'xlim'); xdd = xl(1)-(diff(xl)*.1);
            yl = get(gca,'ylim'); ydd = mean(yl)+[-40 40];
            view(90,90); axis square equal
            xlabel('North >'); ylabel('East >'); zlabel('Down >');
            tmp = get(ovpath_fig,'Position');
            title(['Session ',metadata(i).session]);
            set(ovpath_fig,...
                'DefaultAxesFontSize',5,...
                'DefaultLineLineWidth',2,...
                'PaperUnits','points',...
                'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
                'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
            saveas(ovpath_fig,[pathfigFol 'OverViewScale_' saveL '.pdf'],'pdf');

            % set up for no-vertical plot with stridestats
            kk = 1;
            stridestats = [];
            if figureson ==0
                novert_fig = figure('Position',[0 0 1200 400],'visible','off');
            else
                novert_fig = figure('Position',[0 0 1200 400]);
            end
            cmp = lines(sum(kpp));
            scatter3(posPlot(:,1),posPlot(:,2),posPlot(:,3),'.k')

            flag = 0;
            % more plots in here, decide what plots we want
            for ki = find(kpp)'
                if ki == length(stationaryStart) && kpp(end) == 1
                    ki = ki-1;
                    flag = 1; % set flag so after the loop you can remove the last one
                end
                line(posPlot([stationaryStart(ki)-1 stationaryStart(ki+1)-1],1),...
                    posPlot([stationaryStart(ki)-1 stationaryStart(ki+1)-1],2),...
                    posPlot([stationaryStart(ki)-1 stationaryStart(ki)-1],3),... keep same vertical position
                    'Color',cmp(kk,:),'LineWidth',2)

                stmp = posPlot(stationaryStart(ki+1)-1, :) - posPlot(stationaryStart(ki),:);
                dtmp = (stationaryStart(ki+1)-1 - stationaryStart(ki));

                text(posPlot(stationaryStart(ki),1),...
                    posPlot(stationaryStart(ki),2),...
                    posPlot(stationaryStart(ki),3)+.1,...
                    [num2str(norm(stmp(1:2))) ' ft'],'Color',cmp(kk,:),...
                    'Rotation',-90,'HorizontalAlignment','right')
                text(posPlot(stationaryStart(ki),1),...
                    posPlot(stationaryStart(ki),2),...
                    posPlot(stationaryStart(ki),3)-.1,...
                    [num2str(dtmp*samplePeriod) 's'],'Color',cmp(kk,:),'Rotation',-90);


                stridestats.SL(kk) = norm(stmp); %consider only change in AP-ML plane - changed to 3 planes for sc12 and sc8 analysis
                stridestats.DU(kk) = dtmp*samplePeriod;
                kk=kk+1;
            end

            % more plotting, decide what we want to plot
            zlim(get(gca,'zlim')+[-.3 .3]);
            view(-45,80);
            tmp = get(novert_fig,'Position');
            set(novert_fig,...
                'DefaultAxesFontSize',5,...
                'DefaultLineLineWidth',2,...
                'PaperUnits','points',...
                'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
                'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
            title(['Session ',metadata(i).session]);
            saveas(novert_fig,[pathfigFol 'M2stats_novert_' saveL '.pdf'],'pdf');



            % if flag is true, remove the last SL and DU from stridestats
            if flag == 1
                stridestats.SL = stridestats.SL(1:end-1);
                stridestats.DU = stridestats.DU(1:end-1);
            end

            %If stridestats is empty (as in the case of data with no
            %walking segments), create empty DU and SL fields
            if isempty(stridestats)
                disp(['no walkseg, skipping stride extraction for ' metadata(i).session]);
            else
                % remove short steps and the outliers from stridestats here BEFORE
                % (method from error ellipse)
                % create new version without removing outliers to save later
                stridestats_all = stridestats;
                tmp = zeros(1,length(stridestats.SL));
                tmp(stridestats.SL<.1)=1; % remove small steps
                stridestats.SL = stridestats.SL(~tmp);
                stridestats.DU = stridestats.DU(~tmp);
                eh = error_ellipse(stridestats.SL',stridestats.DU',0);
                tmp2 = eh.outlier;
                stridestats.SL = stridestats.SL(~tmp2); % remove outliers from error ellipse
                stridestats.DU = stridestats.DU(~tmp2);

                % outliers removed version
                stridestats.SLmean = mean(stridestats.SL(stridestats.SL>0.1&stridestats.DU>0.1));
                stridestats.SLsd = std(stridestats.SL(stridestats.SL>0.1&stridestats.DU>0.1));
                stridestats.DUmean = mean(stridestats.DU(stridestats.SL>0.1&stridestats.DU>0.1));
                stridestats.DUsd = std(stridestats.DU(stridestats.SL>0.1&stridestats.DU>0.1));
                stridestats.WSmean = mean(stridestats.SL(stridestats.SL>0.1&stridestats.DU>0.1)./...
                    stridestats.DU(stridestats.SL>0.1&stridestats.DU>0.1));

                % version with all the data
                stridestats_all.SLmean = mean(stridestats_all.SL);
                stridestats_all.SLsd = std(stridestats_all.SL);
                stridestats_all.DUmean = mean(stridestats_all.DU);
                stridestats_all.DUsd = std(stridestats_all.DU);
                stridestats_all.WSmean = mean(stridestats_all.SL./stridestats_all.DU);


                save([stepFol 'M2parameters_' saveL '.mat'],...
                    'stridestats','stridestats_all','posPlot','stationary')

            end
        end
    end
end

end

