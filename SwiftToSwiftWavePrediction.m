clc
clearvars
close all

swiftDirectory = '/Users/mike/Documents/UW/Research/Data/LC-DRI data';
swiftIDs = {'22','23','24','25'};
%swiftDate = '28Mar2017'; swiftTime = '19_02';
swiftDate = '02Apr2017'; swiftDate0 = '02Apr2017';  swiftTime = '23_02';

performSurfaceReconstruction = true;
performPointComparison = false;

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
    minLat(i) = min(imuData(i).GpsPos.lat);
    maxLat(i) = max(imuData(i).GpsPos.lat);
    medLat(i) = median(imuData(i).GpsPos.lat);
    minLon(i) = min(imuData(i).GpsPos.long);
    maxLon(i) = max(imuData(i).GpsPos.long);
    medLon(i) = median(imuData(i).GpsPos.long);
    minTime(i) = min(imuData(i).ShipMotion.time_stamp);
    maxTime(i) = max(imuData(i).ShipMotion.time_stamp);
    time_stamp_step(i) = imuData(i).ShipMotion.time_stamp(end)- imuData(i).ShipMotion.time_stamp(end-1);
    zone = utmzone([minLat(i),minLon(i)]);
    [ellipsoid,estr] = utmgeoid(zone);
    utmstruct = defaultm('utm');
    utmstruct.zone = zone;
    utmstruct.geoid = ellipsoid(1,:);
    utmstruct = defaultm(utmstruct);
    %[utm(i).x,utm(i).y] = mfwdtran(utmstruct,imuData(i).EkfNav.latitude,imuData(i).EkfNav.longitude);
    [utm(i).x,utm(i).y] = mfwdtran(utmstruct,imuData(i).GpsPos.lat,imuData(i).GpsPos.long);
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

t = min(minTime):time_stamp_step:max(maxTime);
nt = length(t);
x_swift = nan(nt,numSwifts);
y_swift = nan(nt,numSwifts);
z_swift = nan(nt,numSwifts);
u_swift = nan(nt,numSwifts);
v_swift = nan(nt,numSwifts);
for i = 1:length(imuData)
%     [tUnique_EkfNav,indUnique_EkfNav,~] = unique(imuData(i).EkfNav.time_stamp);
    [tUnique_GpsPos,indUnique_GpsPos,~] = unique(imuData(i).GpsPos.time_stamp);
    x_swift(:,i) = interp1(tUnique_GpsPos,utm(i).x(indUnique_GpsPos),t,'linear',NaN);
    y_swift(:,i) = interp1(tUnique_GpsPos,utm(i).y(indUnique_GpsPos),t,'linear',NaN);
    [tUnique_GpsVel,indUnique_GpsVel,~] = unique(imuData(i).GpsPos.time_stamp);
    u_swift(:,i) = interp1(tUnique_GpsVel,imuData(i).GpsVel.vel_e(indUnique_GpsVel),t,'linear',NaN);
    v_swift(:,i) = interp1(tUnique_GpsVel,imuData(i).GpsVel.vel_n(indUnique_GpsVel),t,'linear',NaN); 
    [tUnique_ShipMotion,indUnique_ShipMotion,~] = unique(imuData(i).ShipMotion.time_stamp);
    z_swift(:,i) = interp1(tUnique_ShipMotion,imuData(i).ShipMotion.heave(indUnique_ShipMotion),t,'linear',NaN);
    goodInd = t>0.5e8;
    [Hs(i), Tp(i), Dp(i), E(i,:), f(i,:), a1(i,:), b1(i,:), a2(i,:), b2(i,:)] = GPSwaves(u_swift(goodInd,i),v_swift(goodInd,i),z_swift(goodInd,i),5);
    [NS(:,:,i),NE(:,:,i)] = MEM_directionalestimator(a1(i,:),a2(i,:),b1(i,:),b2(i,:),E(i,:),0);
end
f = f(1,:);
theta = 0:2:358;
%%
fig(3) = figure(3); clf(fig(3))
hold on
for i = 1:numSwifts
    plot(f,E(i,:),'-k')
end
hold off
set(gca,'XScale','linear','YScale','linear')
xlim([0 0.6])
%ylim([0 5])
xlabel('f (Hz)')
ylabel('E (m^2/Hz)')

fig(31) = figure(31); clf(fig(31))
pcolor(theta,f,log10(mean(NE,3)));
shading('flat')
set(gca,'XScale','linear','YScale','linear')
xlabel('\theta (deg)')
ylabel('f (Hz)')

%%
fig(4) = figure(4); clf(fig(4));
hold on
for i = 1:numSwifts
    scatter(x_swift(:,i),y_swift(:,i),6,t/1e6)
    %scatter(imuData(i).GpsPos.long,imuData(i).GpsPos.lat,6,imuData(i).GpsPos.time_stamp)
    text(x_swift(end-2,i),y_swift(end-2,i),sprintf('SWIFT %d',i))
end
hold off
xlim([-200 200])
ylim([-200 200])
%xlim(mean(medLon)+[-0.001,0.001])
%ylim(mean(medLat)+[-0.001,0.001])
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
T_delay = -15; %sec
overlap = 1/6;
dt = 0.2; 

% set least squares model values
k = logspace(-3,-1,31)*2*pi;
theta_wavenumber = linspace(-180,180,13)*pi/180;
theta_wavenumber = theta_wavenumber(2:end);
reg_factor = 5e0;

if performPointComparison
    
    for i = 1:numSwifts
        x_array = x_swift(:,setdiff(1:numSwifts,i));
        y_array = y_swift(:,setdiff(1:numSwifts,i));
        z_array = z_swift(:,setdiff(1:numSwifts,i));
        x_target = x_swift(:,i);
        y_target = y_swift(:,i);
        z_target = z_swift(:,i);
        
        [z_pred(:,:,i),z_truth(:,:,i),t_pred(:,:,i)] = runLeastSquaresPrediction_PointComparison(...
            x_target,y_target,z_target,dt,x_array,y_array,z_array,...
            k,theta_wavenumber,reg_factor,T_meas,T_pred,T_delay,overlap);
    end
    
    
    %%
        H_sig = nanmean(Hs);

    fig(6) = figure(6); clf(fig(6));
    fig(7) = figure(7); clf(fig(7));
    for i = 1:numSwifts
        t_plot = t_pred(:,1:round(T_pred*overlap/dt),i)';
        z_pred_plot = z_pred(:,1:round(T_pred*overlap/dt),i)';
        z_truth_plot = z_truth(:,1:round(T_pred*overlap/dt),i)';
        t_plot = t_plot(~isnan(z_pred_plot));
        z_truth_plot = z_truth_plot(~isnan(z_pred_plot));
        z_pred_plot = z_pred_plot(~isnan(z_pred_plot));
        
        
        H_pred = imag(hilbert(z_pred_plot(:)'))';
        phi_pred = atan2(H_pred,z_pred_plot(:));
        A_pred = sqrt(z_pred_plot(:).^2+H_pred.^2);
        H_truth = imag(hilbert(z_truth_plot(:)'))';
        phi_truth = atan2(H_truth,z_pred_plot(:));
        A_truth = sqrt(z_truth_plot(:).^2+H_truth.^2);
        %A_pred = smooth(A_pred,25);
        %A_truth = smooth(A_truth,25);
        
        figure(fig(6))
        subplot(4,1,i)
        hold on
        %plot(t_plot(:),z_pred_plot(:),'-k')
        %plot(t_plot(:),z_truth_plot(:),'-b')
%        plot((0:length(t)-1)*dt,z_swift(:,i),'-')
        plot(t_plot(:),A_pred(:),'-k')
        %plot(t_plot(:),-A_pred(:),'-k')
        plot(t_plot(:),A_truth(:),'-b')
        %plot(t_plot(:),-A_truth(:),'-b')
        
        hold off
        ylim([-H_sig H_sig])
        xlim([100 400])
        legend('Predicted','Measured')
        ylabel('\eta (m)')
        title(sprintf('Swift %d',i))
    xlabel('t (s)')
    
    minT = 100;
    maxT = 450;
        figure(fig(7))
        subplot(2,2,i)
        hold on
        %plot(z_truth_plot(t_plot>minT),z_pred_plot(t_plot>minT),'.')
        %[z_regress,z_regress_int,z_regress_resid,~,regress_stats] = regress(z_truth_plot(t_plot>minT),z_pred_plot(t_plot>minT));
        A_truth_plot = A_truth(t_plot>minT & t_plot<maxT);
        A_pred_plot = A_pred(t_plot>minT & t_plot<maxT);
        plot(A_truth_plot,A_pred_plot,'.')
        [z_regress,z_regress_int,z_regress_resid] = regress(A_pred_plot,A_truth_plot);
        regress_r2 = 1-nansum(z_regress_resid.^2)/nansum((A_truth_plot-nanmean(A_truth_plot)).^2)
        plot([0 H_sig],[0 H_sig],'--k')
        plot([0 H_sig],[0 H_sig].*z_regress,'--b')
        hold off
        xlim([0 H_sig])
        ylim([0 H_sig])
        xlabel('Measured')
        ylabel('Predicted')
        title(sprintf('Swift %d',i))
    end
    
    %%
    fig(9) = figure(9); clf(fig(9));
    i = 1;
        t_plot = t_pred(:,1:round(T_pred*overlap/dt),i)';
        z_pred_plot = z_pred(:,1:round(T_pred*overlap/dt),i)';
        z_truth_plot = z_truth(:,1:round(T_pred*overlap/dt),i)';
        t_plot = t_plot(~isnan(z_pred_plot));
        z_truth_plot = z_truth_plot(~isnan(z_pred_plot));
        z_pred_plot = z_pred_plot(~isnan(z_pred_plot));
        
        subplot(2,2,1:2)
        hold on
        plot(t_plot(:),z_pred_plot(:),'-k')
        plot(t_plot(:),z_truth_plot(:),'-b') 
        hold off
        ylim([-H_sig H_sig])
        xlim([200 300])
        legend('Predicted','Measured')
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
    burst_start = find(t_pred(:,1)>100,1,'first');
    for j = 1:size(z_target_pred,1)
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

