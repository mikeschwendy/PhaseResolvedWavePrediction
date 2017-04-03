clc
clearvars
close all

swiftDirectory = '/Users/mike/Documents/UW/Research/Data/SanNicolasIslandSwift/29March2017_22_01';
resultsDirectory = '/Users/mike/Documents/UW/Research/Data/SanNicolasIslandSwift/29March2017_22_01';

swiftFiles = dir([swiftDirectory '/*.mat']);
numSwifts = length(swiftFiles);
for i = 1:numSwifts
    load([swiftDirectory '/' swiftFiles(i).name]);
    imuData(i) = sbgData;
    minLat(i) = min(imuData(i).EkfNav.latitude);
    maxLat(i) = max(imuData(i).EkfNav.latitude);
    minLon(i) = min(imuData(i).EkfNav.longitude);
    maxLon(i) = max(imuData(i).EkfNav.longitude);
    minTime(i) = min(imuData(i).ShipMotion.time_stamp);
    maxTime(i) = max(imuData(i).ShipMotion.time_stamp);
    dt(i) = imuData(i).ShipMotion.time_stamp(end)- imuData(i).ShipMotion.time_stamp(end-1);
    zone = utmzone([minLat(i),minLon(i)]);
    [ellipsoid,estr] = utmgeoid(zone);
    utmstruct = defaultm('utm');
    utmstruct.zone = zone;
    utmstruct.geoid = ellipsoid(1,:);
    utmstruct = defaultm(utmstruct);
    [utm(i).x,utm(i).y] = mfwdtran(utmstruct,imuData(i).EkfNav.latitude,imuData(i).EkfNav.longitude);
    minX(i) = min(utm(i).x);
    maxX(i) = max(utm(i).x);
    medX(i) = median(utm(i).x);
    minY(i) = min(utm(i).y);
    maxY(i) = max(utm(i).y);
    medY(i) = median(utm(i).y);
end
meanX_tot = mean(medX);
meanY_tot = mean(medY);

fig(1) = figure(1);
subplot(3,1,1)
hold on
for i = 1:numSwifts
    plot(imuData(i).ShipMotion.time_stamp,imuData(i).ShipMotion.heave)
end
hold off
subplot(3,1,2)
hold on
% for i = 1:numSwifts
%     %plot(imuData(i).EkfNav.time_stamp,utm(i).x-meanX_tot)
%     plot(imuData(i).EkfNav.time_stamp,imuData(i).EkfNav.latitude)
% end
%   plot(imuData(1).EkfNav.time_stamp,imuData(1).EkfNav.latitude)
%    plot(imuData(2).EkfNav.time_stamp,imuData(2).EkfNav.latitude,'--')
plot(imuData(3).EkfNav.time_stamp,imuData(3).EkfNav.latitude)
plot(imuData(4).EkfNav.time_stamp,imuData(4).EkfNav.latitude,'--')
hold off
%ylim([-500,500])

subplot(3,1,3)
hold on
% for i = 1:numSwifts
%     %plot(imuData(i).EkfNav.time_stamp,utm(i).y-meanY_tot)
%     plot(imuData(i).EkfNav.time_stamp,imuData(i).EkfNav.longitude)
% end
%    plot(imuData(1).EkfNav.time_stamp,imuData(1).EkfNav.longitude)
%    plot(imuData(2).EkfNav.time_stamp,imuData(2).EkfNav.longitude,'--')
%plot(imuData(4).EkfNav.time_stamp,imuData(4).EkfNav.longitude,'--')
hold off
ylim([0 10])
%ylim([-500,500])

nLat = 101;
nLon = 101;
latVec = linspace(min(minLat),max(maxLat),nLat);
lonVec = linspace(min(minLon),max(maxLon),nLon);
t = minTime:dt:maxTime;
nt = length(t);
x = nan(numSwifts,nt);
y = nan(numSwifts,nt);
z = nan(numSwifts,nt);
for i = 1:length(imuData)
    [tUnique,indUnique,~] = unique(imuData(i).EkfNav.time_stamp);
    x(i,:) = interp1(tUnique,utm(i).x(indUnique),t,'linear',NaN);
    y(i,:) = interp1(tUnique,utm(i).y(indUnique),t,'linear',NaN);
    [tUnique2,indUnique2,~] = unique(imuData(i).ShipMotion.time_stamp);
    z(i,:) = interp1(tUnique2,imuData(i).ShipMotion.heave(indUnique2),t,'linear',NaN);
end

fig(2) = figure(2);
subplot(2,1,1)
hold on
plot(t,x(1,:)-x(2,:))
plot(t,x(3,:)-x(4,:))
plot(t,x(1,:)-x(4,:))
hold off
%ylim([-40 40])
subplot(2,1,2)
hold on
plot(t,y(1,:)-y(2,:))
plot(t,y(3,:)-y(4,:))
plot(t,y(1,:)-y(4,:))
hold off
%ylim([-40 40])

% fig(2) = figure(2);
% for i = 1:nt
%     scatter(x(:,i)-medX,y(:,i),40,z(:,i))
%     xlim([3.2e5,3.23e5])
%     %ylim([3.665e6,3.675e6])
%     axis equal
%     pause(.01);
% end

for i = 1:numSwifts
    goodInd = imuData(i).ShipMotion.time_stamp>0.5e8;
    [S_z(i,:),f_z(i,:)] = pwelch(imuData(i).ShipMotion.heave(goodInd),[],[],512,5);
    goodInd = imuData(i).EkfNav.time_stamp>0.5e8;
    [S_u(i,:),f_u(i,:)] = pwelch(imuData(i).EkfNav.velocity_e(goodInd),[],[],512,5);
    [S_v(i,:),f_v(i,:)] = pwelch(imuData(i).EkfNav.velocity_n(goodInd),[],[],512,5);
end
%%
fig(3) = figure(3); clf(fig(3))
%subplot(3,1,1)
hold on
for i = 1:numSwifts
       plot(f_z(i,:),S_z(i,:),'-k') 
end
%hold off
%set(gca,'XScale','log','YScale','log')
%set(gca,'XScale','linear','YScale','linear')
%subplot(3,1,2)
%hold on
for i = 1:numSwifts
       plot(f_u(i,:),S_u(i,:)./((2*pi*f_u).^2),'-b') 
end
%hold off
%set(gca,'XScale','log','YScale','log')
%subplot(3,1,3)
%hold on
for i = 1:numSwifts
       plot(f_v(i,:),S_v(i,:)./((2*pi*f_v).^2),'-r') 
end
hold off
set(gca,'XScale','log','YScale','log')
%set(gca,'XScale','linear','YScale','linear')
