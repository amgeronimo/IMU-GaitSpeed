%Manually choose segments of continuous walking from the acceleration trace
function [] = choose_walking_segments(numericdata,metadata,printfigure,redcap_path,proc_path,filesep)
 bigOrg = redcap_path.bigOrg;
 rawall = redcap_path.rawall;
 indiv = redcap_path.indiv;
 group = redcap_path.group;

 path_csv = proc_path.csv;
 path_fig = proc_path.fig;
 path_datafig = proc_path.datafig;
 path_pathfig = proc_path.pathfig;
 path_imu = proc_path.imu;
 path_step = proc_path.step;
 path_line = proc_path.line;

    if ispc
        filesep = '\';             % check computer type for file directory 
     elseif ismac                  % separator: / or \
        filesep = '/';
    end
 
%% Extract multiple segments of walking from acceleration and gyroscope 
for i = 1:length(numericdata)
    
    imuFol = [bigOrg,filesep,indiv,filesep,metadata(i).SCname,filesep,metadata(i).SCname,path_imu,filesep];
    figFol = [bigOrg,filesep,indiv,filesep,metadata(i).SCname,filesep,metadata(i).SCname,path_fig,filesep];

    if sum(strcmp(metadata(i).datatype,{'Linear Acceleration','Accelerometer','ACC'}))>0
        if exist([imuFol metadata(i).name '_DATAseg_' metadata(i).date '.mat'],'file')==2
            disp([metadata(i).name '_DATAseg_' metadata(i).date '.mat <strong>exists</strong>'])
            if printfigure == 1
                load(join([imuFol, metadata(i).SCname '_DATAseg_' metadata(i).date '.mat']))
                gg = figure('Visible','off','Position',[0 0 1200 400]); hold on;
                tvec = numericdata(i).time_interp_shift;
                plot(tvec, numericdata(i).pts_interp_shift);
                yll = get(gca,'ylim');
                cmp = parula(size(walkseg,1));
                for jx = 1:size(walkseg,1)
                    ttmp = tvec([walkseg(jx,1) walkseg(jx,2)]);

                    patch([ttmp(1) ttmp(2) ttmp(2) ttmp(1) ttmp(1)],...
                        [yll(1) yll(1) yll(2) yll(2) yll(1)],cmp(jx,:),...
                        'HandleVisibility','off','EdgeAlpha',0,'FaceAlpha',.1);
                end
                xlabel('Time (s)');
                tmp = get(gg,'Position');
                set(gg,...
                    'DefaultAxesFontSize',5,...
                    'DefaultLineLineWidth',2,...
                    'PaperUnits','points',...
                    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
                    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
                saveas(gg,join([figFol, metadata(i).SCname, path_pathfig, metadata(i).SCname '_Routesegs_' metadata(i).date '.pdf']),'pdf');
                close(gg);
            end
        else

            disp(['Processing Patient ' metadata(i).SCname])
            gg=figure('Position',[0 100 1500 600],'Name',metadata(i).SCname); hold on;
            ixx = 1;
            cmp = parula(2);
            for ix = [i i+1]
                plot(numericdata(ix).time_interp_shift, numericdata(ix).pts_interp_shift,'color',cmp(ixx,:));
                ixx = ixx+1;
            end
            yll = get(gca,'ylim');

            xlabel('Time (s)');
            cmp = parula(8);
            pind = 1;
            kk=1;
            walkseg = [];
            while 1
                title('Click two points to select a WALKING bounding window. If none, press enter.');
                [x,~] = ginput(2);
                if ~isempty(x)
                    [~,xi] = min(abs(numericdata(ix).time_interp_shift-x(1)));
                    [~,xe] = min(abs(numericdata(ix).time_interp_shift-x(2)));
                    walkseg(kk,:) = [xi xe];
                    h=patch([x(1) x(2) x(2) x(1) x(1)],[yll(1) yll(1) yll(2) yll(2) yll(1)],cmp(pind,:),'FaceAlpha',.2,'EdgeAlpha',0);
                end
                
                title('Click the mouse to define another WALKING window. If finished, press enter');
                sl = waitforbuttonpress;
                kk=kk+1;
                
                if pind == 4
                    pind = 1;
                else
                    pind = pind+1;
                end
                
                if sl == 1
                    break
                end
            end
            
            save([imuFol,metadata(i).name,'_DATAseg_',metadata(i).date,'.mat'],'walkseg');                      
        end
    end
end
