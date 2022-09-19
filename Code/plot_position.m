function [] = plot_position(NM,acct,samplePeriod,stationary,walkseg)
% acc is in m/s^2

% Integrate acceleration to yield velocity
vel = zeros(size(acct));
for t = 2:length(vel)
    vel(t,:) = vel(t-1,:) + acct(t,:) * samplePeriod;
    if(stationary(t) == 1)
        vel(t,:) = [0 0 0];     % force zero velocity when foot stationary
    end
end
% missing cumtrapz

%         % Compute integral drift during non-stationary periods
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
    for i = 1:numel(stationaryEnd)
        driftRate = vel(stationaryEnd(i)-1, :) / (stationaryEnd(i) - stationaryStart(i));
        enum = 1:(stationaryEnd(i) - stationaryStart(i));
        drift = [enum'*driftRate(1) enum'*driftRate(2) enum'*driftRate(3)];
        velDrift(stationaryStart(i):stationaryEnd(i)-1, :) = drift;
    end
else
    for i = 1:numel(stationaryEnd)
        driftRate = vel(stationaryEnd(si)-1, :) / (stationaryEnd(si) - stationaryStart(si));
        enum = 1:(stationaryEnd(si) - stationaryStart(si));
        drift = [enum'*driftRate(1) enum'*driftRate(2) enum'*driftRate(3)];
        velDrift(stationaryStart(si):stationaryEnd(si)-1, :) = drift;
    end
end

% Remove integral drift
vel = vel - velDrift;

% -------------------------------------------------------------------------
% Compute translational position

% Integrate velocity to yield position
pos = zeros(size(vel));
for t = 2:length(pos)
    pos(t,:) = pos(t-1,:) + vel(t,:) * samplePeriod;    % integrate velocity to yield position
end


% if isempty(stationaryStart)
%     return;
% end
% posPlot = pos(stationaryStart(1):stationaryStart(end),:); 
% ws = walkseg-stationaryEnd(1)+1;


posPlot = pos(min(walkseg(:)):max(walkseg(:)),:);
ws = walkseg - min(walkseg(:)) + 1;

cmp = winter(size(posPlot,1));

subplot(221);
scatter3(posPlot(:,1),posPlot(:,2),posPlot(:,3),10,cmp,'Marker','.')
view(90,90); axis square equal
xlabel('Y >'); ylabel('X >'); zlabel('Z >');

subplot(222);
scatter3(posPlot(:,1),posPlot(:,2),posPlot(:,3),10,cmp,'Marker','.')
view(45,45);
xlabel('Y >'); ylabel('X >'); zlabel('Z >');

subplot(223);
scatter3(posPlot(:,1),posPlot(:,2),posPlot(:,3),10,cmp,'Marker','.')
view(90,0); axis square equal
xlabel('Y >'); ylabel('X >'); zlabel('Z >');

subplot(224);
scatter3(posPlot(:,1),posPlot(:,2),posPlot(:,3),10,cmp,'Marker','.')
view(0,0); axis square equal
xlabel('Y >'); ylabel('X >'); zlabel('Z >');
title('All units in meters')


if isempty(ws)
    return;
end
% % check if walking segment ends after the last stationaryStart - if it
% % does, it will cause an error when highlighting the WS in red
% if ws(end,end) > stationaryStart(end)-stationaryStart(1)
%     % make new WS variable to not overwrite walkseg
%     ws(end,end) = stationaryStart(end)-stationaryStart(1); % shift so the end of the ws is right before data ends
% end
% % check if walking segment starts before the first stationaryStart
% if ws(1,1) < stationaryStart(1)
%     ws(1,1) = stationaryStart(1);
% end
% % check if any of the walking segments end after the last stationaryStart
% for w = 1:size(ws,1)
%     % check end of ws
%     if ws(w,2) > stationaryStart(end)-stationaryStart(1)
%         ws(w,2) = stationaryStart(end)-stationaryStart(1);
%     else
%         % check beginning of ws
%         if ws(w,1) < stationaryStart(1)
%             ws(1,1) = stationaryStart(1);
%         end
%     end
% end

% add walk seg in red
for i = 1:size(ws)
    subplot(221)
    hold on
    scatter3(posPlot(ws(i,1):ws(i,2),1),posPlot(ws(i,1):ws(i,2),2),posPlot(ws(i,1):ws(i,2),3),10,'r','Marker','.')
    hold off
    subplot(222)
    hold on
    scatter3(posPlot(ws(i,1):ws(i,2),1),posPlot(ws(i,1):ws(i,2),2),posPlot(ws(i,1):ws(i,2),3),10,'r','Marker','.')
    hold off
    subplot(223)
    hold on
    scatter3(posPlot(ws(i,1):ws(i,2),1),posPlot(ws(i,1):ws(i,2),2),posPlot(ws(i,1):ws(i,2),3),10,'r','Marker','.')
    hold off
    subplot(224)
    hold on
    scatter3(posPlot(ws(i,1):ws(i,2),1),posPlot(ws(i,1):ws(i,2),2),posPlot(ws(i,1):ws(i,2),3),10,'r','Marker','.')
    hold off
end
