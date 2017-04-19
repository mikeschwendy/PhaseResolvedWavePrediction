clc
clearvars
close all

swiftDirectory = '/Users/mike/Documents/UW/Research/Data/LC-DRI data';
swiftIDs = {'22','23','24','25'};
%swiftDate = '28Mar2017'; swiftTime = '19_02';
swiftDate = '02Apr2017'; swiftDate0 = '02Apr2017';  swiftTime = '23_04';

performSurfaceReconstruction = false;
performPointComparison = true;

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
x_swift = nan(nt,numSwifts);
y_swift = nan(nt,numSwifts);
z_swift = nan(nt,numSwifts);
for i = 1:length(imuData)
    [tUnique,indUnique,~] = unique(imuData(i).EkfNav.time_stamp);
    x_swift(:,i) = interp1(tUnique,utm(i).x(indUnique),t,'linear',NaN);
    y_swift(:,i) = interp1(tUnique,utm(i).y(indUnique),t,'linear',NaN);
    [tUnique2,indUnique2,~] = unique(imuData(i).ShipMotion.time_stamp);
    z_swift(:,i) = interp1(tUnique2,imuData(i).ShipMotion.heave(indUnique2),t,'linear',NaN);
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
set(gca,'XScale','linear','YScale','linear')
xlim([0 0.6])
xlabel('f (Hz)')
ylabel('E (m^2/Hz)')

%%
fig(4) = figure(4); clf(fig(4));
hold on
for i = 1:numSwifts
    scatter(utm(i).x,utm(i).y,6,imuData(i).EkfNav.time_stamp/1e6)
end
hold off
xlim([-200 200])
ylim([-200 200])
xlabel('x (m)')
ylabel('y (m)')
cbar = colorbar;
ylabel(cbar,'t (s)')

%%
fig(5) = figure(5); clf(fig(5));
for i = 1:numSwifts
    subplot(4,1,i)
    plot(t/1e6,z_swift(:,i),'k')
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
k = logspace(-3,-1,31)*2*pi;
theta_wavenumber = linspace(-180,180,19)*pi/180;
theta_wavenumber = theta_wavenumber(2:end);
reg_factor = 1e2;

if performPointComparison
    
    for i = 1:numSwifts
        x_array = x_swift(:,setdiff(1:numSwifts,i));
        y_array = y_swift(:,setdiff(1:numSwifts,i));
        z_array = z_swift(:,setdiff(1:numSwifts,i));
        x_target = x_swift(:,i);
        y_target = y_swift(:,i);
        z_target = z_swift(:,i);
        
        [z_pred(:,:,i),z_truth(:,:,i),t_pred(:,:,i)] = runLeastSquaresPrediction_PointComparison(...
            x_target,y_target,z_target,0.2,x_array,y_array,z_array,...
            k,theta_wavenumber,reg_factor,T_meas,T_pred,T_delay,overlap);
    end
    
    H_sig = 1;
    
    fig(6) = figure(6); clf(fig(6));
    for i = 1:numSwifts
        subplot(4,1,i)
        hold on
        plot(t_pred(:,:,i)',z_pred(:,:,i)','-k')
        plot((0:length(t)-1)*0.2,z_swift(:,i),'-')
        hold off
        ylim([-H_sig H_sig])
        xlim([100 400])
        legend('Predicted','Measured')
        ylabel('\eta (m)')
        title(sprintf('Swift %d',i))
    end
    xlabel('t (s)')
    
    minT = 100;
    fig(7) = figure(7); clf(fig(7));
    for i = 1:numSwifts
        t_pred_i = t_pred(:,:,i);
        z_truth_i = z_truth(:,:,i);
        z_pred_i = z_pred(:,:,i);
        subplot(2,2,i)
        hold on
        plot(z_truth_i(t_pred_i>minT),z_pred_i(t_pred_i>minT),'.')
        [z_regress,z_regress_int,z_regress_resid,~,regress_stats] = regress(z_truth_i(:),z_pred_i(:));
        regress_r2 = regress_stats(1);
        plot([-H_sig H_sig],[-H_sig H_sig],'--k')
        plot([-H_sig H_sig],[-H_sig H_sig].*z_regress,'--b')
        hold off
        xlim([-H_sig H_sig])
        ylim([-H_sig H_sig])
        xlabel('Measured')
        ylabel('Predicted')
        title(sprintf('Swift %d',i))
    end
end

%%
if performSurfaceReconstruction
    [x_target,y_target] = meshgrid(linspace(-100,100,20),linspace(-100,100,21));
    [ny_target,nx_target] = size(x_target);
    [z_target_pred,t_pred] = runLeastSquaresPrediction_SurfaceReconstruction(...
        x_target,y_target,x_swift,y_swift,z_swift,0.2,...
        k,theta_wavenumber,reg_factor,T_meas,T_pred,T_delay,overlap);
    [num_bursts,Nt_pred,~] = size(z_target_pred);
    z_target_pred = reshape(z_target_pred,[num_bursts,Nt_pred,ny_target,nx_target]);
    
    
    fig(8) = figure(8); clf(fig(8));
    k = 0;
    for j = 4:size(z_target_pred,1)
        for i = 1:size(z_target_pred,2)
            clf(fig(8));
            subplot(4,1,1:3)
            pcolor(x_target,y_target,squeeze(z_target_pred(j,i,:,:)))
            shading('flat')
            t_i = t_pred(j,i,1);
            t_ind = round(t_i/0.2);
            hold on
            H = scatter(x_swift(t_ind,:),y_swift(t_ind,:),100,z_swift(t_ind,:),'filled');
            H.MarkerEdgeColor = 'k';
            hold off
            set(gca,'CLim',[-1 1],'XLim',[-100 100],'YLim',[-100 100])
            subplot(4,1,4)
            plot((0:nt-1)'*0.2*ones(1,numSwifts),z_swift(:,:))
            hold on
            plot(t_i*ones(numSwifts,1),z_swift(t_ind,:),'o')
            hold off
            set(gca,'YLim',[-1 1],'XLim',[100 200])
            %print('-djpeg',['/Users/mike/Documents/UW/Research/Results/LC_DRI_Results/SurfaceReconstruction/Frame_' sprintf('%03d',k) '.jpg'])
            k = k+1;
            pause(.1)
        end
    end
end

