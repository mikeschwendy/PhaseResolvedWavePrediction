%clc
clearvars
close all

swiftDirectory = '/Users/mike/Documents/UW/Research/Data/LC-DRI data';
swiftIDs = {'22','23','24','25'};
%swiftDate = '28Mar2017'; swiftTime = '19_02';
swiftDate = '02Apr2017'; swiftDate0 = '02Apr2017';  swiftTime = '23_02';



swiftFiles = dir([swiftDirectory '/*.mat']);
numSwifts = length(swiftIDs);
for i = 1:numSwifts
    dirI = [swiftDirectory '/SWIFT' char(swiftIDs(i)) '_' swiftDate0 '/SBG/Raw/' datestr(datenum(swiftDate,'ddmmmyyyy'),'yyyymmdd')];
    fileI = [dirI '/SWIFT' char(swiftIDs(i)) '_SBG_' swiftDate '_' swiftTime '.mat'];
    if ~exist(fileI,'file')
        error('SWIFT file not found')
    end
    load(fileI)
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
for i = 1:numSwifts
    utm(i).x = utm(i).x - meanX_tot;
    utm(i).y = utm(i).y - meanY_tot;
end

t = min(minTime):dt:max(maxTime);
nt = length(t);
x_swifts = nan(nt,numSwifts);
y_swifts = nan(nt,numSwifts);
z_swifts = nan(nt,numSwifts);
for i = 1:length(imuData)
    [tUnique,indUnique,~] = unique(imuData(i).EkfNav.time_stamp);
    x_swifts(:,i) = interp1(tUnique,utm(i).x(indUnique),t,'linear',NaN);
    y_swifts(:,i) = interp1(tUnique,utm(i).y(indUnique),t,'linear',NaN);
    [tUnique2,indUnique2,~] = unique(imuData(i).ShipMotion.time_stamp);
    z_swifts(:,i) = interp1(tUnique2,imuData(i).ShipMotion.heave(indUnique2),t,'linear',NaN);
end

for i = 1:numSwifts
    goodInd = imuData(i).ShipMotion.time_stamp>0.5e8;
    [S_z(i,:),f_z(i,:)] = pwelch(imuData(i).ShipMotion.heave(goodInd),[],[],512,5);
    goodInd = imuData(i).EkfNav.time_stamp>0.5e8;
    [S_u(i,:),f_u(i,:)] = pwelch(imuData(i).EkfNav.velocity_e(goodInd),[],[],512,5);
    [S_v(i,:),f_v(i,:)] = pwelch(imuData(i).EkfNav.velocity_n(goodInd),[],[],512,5);
end
%%
fig(3) = figure(3); clf(fig(3))
hold on
for i = 1:numSwifts
    plot(f_z(i,:),S_z(i,:),'-k')
end
% for i = 1:numSwifts
%     plot(f_u(i,:),S_u(i,:)./((2*pi*f_u).^2),'-b')
% end
% for i = 1:numSwifts
%     plot(f_v(i,:),S_v(i,:)./((2*pi*f_v).^2),'-r')
% end
hold off
set(gca,'XScale','linear','YScale','linear','xlim',[0 0.5],'ylim',[0 10])

%%
fig(4) = figure(4); clf(fig(4));
hold on
for i = 1:numSwifts
    scatter(utm(i).x,utm(i).y,6,imuData(i).EkfNav.time_stamp)
end
hold off
xlim([-100 100])
ylim([-100 100])
xlabel('x [km]')
ylabel('y [km]')

%%
fig(5) = figure(5); clf(fig(5));
for i = 1:numSwifts
    subplot(4,1,i)
    plot(t/1e6,z_swifts(:,i),'k')
    xlim([60 180])
    ylim([-2 2])
end
%%
% set time bursts duration/delay etc
T_meas = 30; %sec
T_pred = 30; %sec
T_delay = -30; %sec
overlap = 1;

% set least squares model values
k = logspace(-3,-1,41)*2*pi;
theta_wavenumber = linspace(-180,180,19)*pi/180;
theta_wavenumber = theta_wavenumber(2:end);
reg_factor = 1e3;

% x_array = x_swift(:,setdiff(1:numSwift,ind_target));
% y_array = y_swift(:,setdiff(1:numSwift,ind_target));
% z_array = z_swift(:,setdiff(1:numSwift,ind_target));
% x_target = x_swift(:,ind_target);
% y_target = y_swift(:,ind_target);
% z_target = nan(size(x_target));

[x_target,y_target] = meshgrid(linspace(-100,100,20),linspace(-100,100,21));
[ny_target,nx_target] = size(x_target);
[z_target_pred,t_pred] = runLeastSquaresPrediction_Swifts(...
    x_target,y_target,x_swifts,y_swifts,z_swifts,0.2,...
    k,theta_wavenumber,reg_factor,T_meas,T_pred,T_delay,overlap);
[num_bursts,Nt_pred,~] = size(z_target_pred);
z_target_pred = reshape(z_target_pred,[num_bursts,Nt_pred,ny_target,nx_target]);

%%
if true
fig(8) = figure(8); clf(fig(8));
for j = 2:size(z_target_pred,1)
    for i = 1:size(z_target_pred,2)
        clf(fig(8));
        subplot(4,1,1:3)
        pcolor(x_target,y_target,squeeze(z_target_pred(j,i,:,:)))
        shading('flat')
        t_i = t_pred(j,i,1);
        t_ind = round(t_i/0.2);
        hold on
        H = scatter(x_swifts(t_ind,:),y_swifts(t_ind,:),100,z_swifts(t_ind,:),'filled');
        H.MarkerEdgeColor = 'k';
        hold off
        set(gca,'CLim',[-1 1],'XLim',[-100 100],'YLim',[-100 100])
        subplot(4,1,4)
        plot((0:nt-1)'*0.2*ones(1,numSwifts),z_swifts(:,:))
        hold on
        plot(t_i*ones(numSwifts,1),z_swifts(t_ind,:),'o')
        hold off
        set(gca,'YLim',[-1 1],'XLim',[100 200])
        pause(.1)
    end
end
end

if false
    fig(6) = figure(6); clf(fig(6));
    for i = 1:numSwifts
        subplot(4,1,i)
        hold on
        plot(t_pred(:,:,i)',z_pred(:,:,i)','-k')
        plot((0:length(t)-1)*0.2,z_swifts(:,i),'-')
        hold off
        ylim([-2 2])
    end
    
    fig(7) = figure(7); clf(fig(7));
    H_sig = 2;
    for i = 1:numSwifts
        z_truth_i = z_truth(:,:,i);
        z_pred_i = z_pred(:,:,i);
        subplot(2,2,i)
        hold on
        plot(z_truth_i(:),z_pred_i(:),'.')
        [z_regress,z_regress_int,z_regress_resid,~,regress_stats] = regress(z_truth_i(:),z_pred_i(:));
        regress_r2 = regress_stats(1);
        plot([-H_sig H_sig],[-H_sig H_sig],'--k')
        plot([-H_sig H_sig],[-H_sig H_sig].*z_regress,'--b')
        hold off
        xlim([-H_sig H_sig])
        ylim([-H_sig H_sig])
    end
end

