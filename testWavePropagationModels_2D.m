clc
clearvars
close all

% Load waverider data from Station Papa Jan. 2015
loadSimulatedData = true;
plotSurface = false;
plotFrequencySpectra = true;
plotRenderedFrames = false;
resultsDir = '/Users/mike/Documents/UW/Research/Results/SurfaceSimulations_2D';
% choose data source as SWIFT, CDIP, or Waverider
dataSource = 'CDIP';

switch dataSource
    case 'Waverider'
        % Load Waverider
        spectraDir = '/Users/mike/Documents/UW/Research/Data/StereoTests';
        allWaverider = load([spectraDir '/Waverider/allDirectionalSpectra.mat']);
        
        % Enter start time
        startTime = datenum(2015,1,5,20,30,0); %datenum(2014,12,28,19,29,0); %
        mkdir([resultsDir '/' datestr(startTime,'ddmmmyyyy')])
        plotDir = [resultsDir '/' datestr(startTime,'ddmmmyyyy/HHMMUTC')];
        plotDir = [plotDir '/LongTimeSmallDomain'];
        %    plotDir = [plotDir '/ShortTimeLargeDomain'];
        mkdir(plotDir)
        
        % Find data timestamped within 30 minutes of start time
        waveriderInds = find(abs(allWaverider.time - startTime) < datenum(0,0,0,0,30,0));
        
        % Average waverider directional spectra
        EfthetaWaverider = nanmean(allWaverider.Sftheta(:,:,waveriderInds),3);
        EfthetaWaverider = pi/180*EfthetaWaverider; % convert to m^2/(Hz*deg)
        thetaWaverider = allWaverider.theta;
        fWaverider = allWaverider.f;
        
        % Transform into WAFO spec struct format
        S_Measured = struct();
        S_Measured.date = startTime;
        S_Measured.type = 'dir';
        S_Measured.S = EfthetaWaverider'*180/pi;
        %S_Waverider.S = zeros(length(thetaWaverider),length(fWaverider));
        %S_Waverider.S((length(thetaWaverider)+1)/2,:) = sum(pi/180*4*EfthetaWaverider,2)'*180/pi;
        %S_Waverider.S(1,:) = sum(pi/180*4*EfthetaWaverider,2)'*180/pi;
        S_Measured.f = fWaverider;
        S_Measured.theta = thetaWaverider'*pi/180;
        S_Measured.phi = 0;
        
    case 'SWIFT'
        % Load SWIFT
        swiftDir = '/Users/mike/Documents/UW/Research/Results/PapaSwiftVisualization';
        swiftInd = '01052';
        load([swiftDir '/SWIFTdata_' swiftInd '.mat'],'SWIFT1','SWIFT2');
        plotDir = [resultsDir '/Papa_SWIFT/SWIFT_' swiftInd];
        mkdir(plotDir)
        [EfthetaSwift, theta, freq] = SWIFTdirectionalspectra(SWIFT1, false);
        theta(theta > 180) = theta(theta > 180) - 360;
        [theta, dsort] = sort(theta);
        EfthetaSwift = EfthetaSwift(:,dsort);
        indNan = sum(isnan(EfthetaSwift),2)>0;
        
        % Transform into WAFO spec struct format
        S_Measured = struct();
        S_Measured.date = SWIFT1.time;
        S_Measured.type = 'dir';
        S_Measured.S = EfthetaSwift(~indNan,:)'*180/pi;
        S_Measured.f = freq(~indNan);
        S_Measured.theta = theta'*pi/180;
        S_Measured.phi = 0;
        
    case 'CDIP'
        % load from CDIP buoy
        cdipDir = '/Users/mike/Dropbox/SanNicolasIsland/CDIP 138 Files';
        cdipFile = [cdipDir '/March_2016.dat'];
        % date/time string MUST be of the form 'YYYYMMDDHHmm'
        %        cdipTimestring = '201603010121';
        %        cdipTimestring = '201603020251';
        %        cdipTimestring = '201603030121';
        %        cdipTimestring = '201603040121';
        %cdipTimestring = '201603040151';
        cdipTimestring = '201603041521';
        [tspc,tsys] = read_cdip_buoy(cdipFile,cdipTimestring);
        [~, Etheta] = MEM_directionalestimator(tspc.a1,tspc.a2,tspc.b1,tspc.b2,tspc.en,0);
        dtheta = 2;
        theta = -(-180:dtheta:179);  % start with cartesion (a1 is positive east velocities, b1 is positive north)
        plotDir = [resultsDir '/SanNicolasIsland_CDIP/' cdipTimestring];
        mkdir(plotDir)
        % change direction from ccw from east to cw from north
        theta = -theta + 90;
        % bring into -180:180 range (necessary for WAFO)
        theta(theta > 180) = theta(theta > 180) - 360;
        [theta, dsort] = sort(theta);
        Etheta = Etheta(:,dsort);
        
        % Transform into WAFO spec struct format
        S_Measured = struct();
        S_Measured.date = datenum(cdipTimestring,'yyyymmddHHMM');
        S_Measured.type = 'dir';
        % be sure directions are in radians and energy is similarly scaled
        S_Measured.S = Etheta'*180/pi;
        S_Measured.f = tspc.fr;
        S_Measured.theta = theta'*pi/180;
        S_Measured.phi = 0;
end
%% Plot directional wave spectrum
% mean and peak direction
dir_mean = cell2mat(dspec2char(S_Measured,'Mdir'))*180/pi;
dir_peak = cell2mat(dspec2char(S_Measured,'Wdir'))*180/pi;
% frequency-dependent mean direction and spread
dir_mean_f = cell2mat(dspec2char(S_Measured,'FMdir'))*180/pi;
dir_spread_f = cell2mat(dspec2char(S_Measured,'FSpr'))*180/pi;
% mean and peak frequency
f_mean = 1./spec2char(S_Measured,'Tm01');
f_peak = 1./spec2char(S_Measured,'Tp');
% significant wave height
H_sig = spec2char(S_Measured,'Hm0');
% mean and peak wavenumber and phase speed from linear dispersion
k_mean = (2*pi*f_mean).^2/9.8;
k_peak = (2*pi*f_peak).^2/9.8;
c_mean = 2*pi*f_mean/k_mean;
c_peak = 2*pi*f_peak/k_peak;
% optional - rotate spectra
thetaRot = S_Measured.theta - dir_peak*pi/180;
thetaRot(thetaRot > pi) = thetaRot(thetaRot > pi) - 2*pi;
thetaRot(thetaRot < -pi) = thetaRot(thetaRot < -pi) + 2*pi;
[thetaRot, dsort] = sort(thetaRot);
S_Measured.S = S_Measured.S(dsort,:);
dir_mean = dir_mean - dir_peak;
dir_peak = 0;

%% Plot directional spectra
fig(6) = figure(6); clf(fig(6));
fig(6).Position = [1 1 5 5];
fig(6).Position = fig(6).PaperPosition;
pcolor(S_Measured.theta*180/pi,S_Measured.f,log10(S_Measured.S'))
%pcolor(S_Measured.theta*180/pi,S_Measured.f,S_Measured.S')
shading('flat')
hold on
plot(dir_peak,f_peak,'xk')
plot(dir_mean,f_mean,'xr')
hold off
xlabel('\theta (degrees)')
ylabel('f (Hz)')
print(gcf,[plotDir '/DirectionalSpectra.png'],'-dpng')

% For reference...
wMax = max(2*pi*S_Measured.f);
wMin = min(2*pi*S_Measured.f);
kMax = w2k(wMax);
kMin = w2k(wMin);
lMin = 2*pi/kMax;  % dx should equal approx. lMin/4
lMax = 2*pi/kMin;
tMin = 2*pi/wMax;  % dt should equal approx. tMin/4
tMax = 2*pi/wMin;

%% Perform Simulation
if ~loadSimulatedData
    % Define simulation domain
    dx = 4;
    dt = 1;
    L = 500;
    T = 300;
    Nx = round(L/dx);
    Ny = Nx+1;
    Nt = round(T/dt);
    x = linspace(0,L-dx,Nx); % m
    y = linspace(0,L,Ny); % m
    t = linspace(0,T-dt,Nt); % s
    dy = dx;
    % Call WAFO codes to simulate surface trace
    % if fftdim = 2, 2D ifft at each time.  if fftdim = 1, 1D ifft at each space
    fftdim = 1;
    [Y,~] = seasim(S_Measured,Nx,Ny,Nt,dx,dy,dt,fftdim,0);
    save([plotDir '/SimulationData_2D.mat'],'Y','dx','dt','L','T','Nx','Nt','x','t','Ny','dy','y','S_Measured')
else
    load([plotDir '/SimulationData_2D.mat'],'Y','dx','dt','L','T','Nx','Nt','x','t','Ny','dy','y','S_Measured')
end

%% Check Frequency Spectra
if plotFrequencySpectra
    binAvg = 1;
    windowLengthSeconds = 120;
    windowLengthSim = 2^nextpow2(windowLengthSeconds/dt);
    fSim = 0:(binAvg/(dt*windowLengthSim)):(1/(2*dt));
    nfSim = length(fSim);
    xSubSim = round(linspace(1,Nx,min(40,Nx)));
    ySubSim = round(linspace(1,Ny,min(40,Ny)));
    SfSim = nan(length(xSubSim),length(ySubSim),nfSim);
    for i = 1:length(xSubSim)
        for j = 1:length(ySubSim)
            % 2*pwelch because giving frequency vector requires two-sided
            SfSim(i,j,:) = 2*pwelch(squeeze(Y.Z(ySubSim(j),xSubSim(i),:)),windowLengthSim,0.5*windowLengthSim,fSim,1/dt,'psd');
        end
    end
    Sf_Measured = spec2spec(S_Measured,'freq');
    
    % Plot Frequency Spectra
    SfSim_avg = squeeze(nanmean(nanmean(SfSim,1),2));
    dtheta = S_Measured.theta(2)-S_Measured.theta(1);
    fig(2) = figure(2); clf(fig(2));
    subplot(2,1,1)
    %plot(Sf_Waverider.f,Sf_Waverider.S,'-k','linewidth',2)
    hold on
    plot(Sf_Measured.f,nansum(S_Measured.S*dtheta),'-k','linewidth',2)
    plot(fSim,SfSim_avg,'-k')
    hold off
    %set(gca,'YScale','log','XScale','log','xlim',[0.05,0.6],'ylim',[1e-4 1e3])
    set(gca,'YScale','linear','XScale','linear','xlim',[0,0.3]);%,'ylim',[0 30])
    legend('Original Waverider Spectrum','Simulation','Stereo','SWIFT',...
        'location','northeast')
    xlabel('f (Hz)')
    ylabel('E(f) (m^2/Hz)')
    subplot(2,1,2)
    plot(Sf_Measured.f,dir_mean_f,'-k','linewidth',2)
    hold on
    plot(Sf_Measured.f,dir_mean_f+dir_spread_f,'-k','linewidth',1)
    plot(Sf_Measured.f,dir_mean_f-dir_spread_f,'-k','linewidth',1)
    hold off
    set(gca,'xlim',[0,0.3],'ylim',[-270,270],'ytick',[-180 -90 0 90 180])
    xlabel('f (Hz)')
    ylabel('\theta (deg)')
    legend('Mean Direction','Spread','location','northeast')
    print(fig(2),[plotDir '/CheckFrequencySpectra.png'],'-dpng','-r300')
    
end

%% Setup buoy array - Circle
if false
    x_target_ind = round(Nx/2);
    y_target_ind = round(Ny/2);
    x_target = x(x_target_ind);
    y_target = y(y_target_ind);
    n_theta = 18;
    theta_spread = 360; %120; % deg
    dir_principal = dir_peak;%dir_mean;%
    theta = linspace(dir_principal-theta_spread/2,dir_principal+theta_spread/2,n_theta)';
    %theta = [dir_principal];
    circle_radius = [50, 100, 150]; % m
    n_radius = length(circle_radius);
    
    circle_radius_vec = repmat(circle_radius(:),[1,n_theta]);
    theta_vec = repmat(theta(:)',[n_radius,1]);
    n_points = n_theta*n_radius;
    x_sentinel = x_target + round(circle_radius_vec(:).*cos(theta_vec(:)*pi/180)/dx)*dx;
    y_sentinel = y_target + round(circle_radius_vec(:).*sin(theta_vec(:)*pi/180)/dy)*dy;
    x_sentinel_ind = nan(n_points,1);
    y_sentinel_ind = nan(n_points,1);
    z_sentinel = nan(n_points,Nt);
    for i = 1:n_points
        x_sentinel_ind(i) = find(Y.x == x_sentinel(i));
        y_sentinel_ind(i) = find(Y.y == y_sentinel(i));
        z_sentinel(i,:) = squeeze(Y.Z(y_sentinel_ind(i),x_sentinel_ind(i),:));
    end
end
%% Setup buoy array - box
if true
    dir_principal = dir_peak;
    box_side_length = 200;
    x_center_ind = round(Nx/2);
    y_center_ind = round(Ny/2);
    x_center = x(x_center_ind);
    y_center = y(y_center_ind);
    x_sentinel = x_center + sqrt(2)*box_side_length/2*[1 0 -1 0 0.5 0.5 0.5 0];
    y_sentinel = y_center + sqrt(2)*box_side_length/2*[0 1 0 -1 0 0.5 -0.5 0];
    n_points = length(x_sentinel);
    R = [cos(dir_principal*pi/180), -sin(dir_principal*pi/180); sin(dir_principal*pi/180), cos(dir_principal*pi/180)];
    xy_sentinel_rot = R*[(x_sentinel-x_center);(y_sentinel-y_center)];
    x_sentinel = round((x_center+xy_sentinel_rot(1,:))/dx)*dx;
    y_sentinel = round((y_center+xy_sentinel_rot(2,:))/dy)*dy;
    
    x_sentinel_ind = nan(n_points,1);
    y_sentinel_ind = nan(n_points,1);
    z_sentinel = nan(n_points,Nt);
    for i = 1:n_points
        x_sentinel_ind(i) = find(Y.x == x_sentinel(i));
        y_sentinel_ind(i) = find(Y.y == y_sentinel(i));
        z_sentinel(i,:) = squeeze(Y.Z(y_sentinel_ind(i),x_sentinel_ind(i),:));
    end
    ind_target = 3;
    ind_target = 1:n_points == ind_target;
    x_target = x_sentinel(ind_target);
    y_target = y_sentinel(ind_target);
    x_target_ind = x_sentinel_ind(ind_target);
    y_target_ind = y_sentinel_ind(ind_target);
    x_sentinel = x_sentinel(~ind_target);
    y_sentinel = y_sentinel(~ind_target);
    z_sentinel = z_sentinel(~ind_target,:);
    x_sentinel_ind = x_sentinel_ind(~ind_target);
    y_sentinel_ind = y_sentinel_ind(~ind_target);
    n_points = n_points-1;
end

%% Setup buoy array - lines


%% Least squares calculation

% set time bursts duration/delay etc
T_meas = 30; %sec
Nt_meas = round(T_meas/dt);
T_pred = 90; %sec
Nt_pred = round(T_pred/dt);
T_delay = 0; %sec
Nt_delay = round(T_delay/dt);
overlap = 1;
t_meas_start_ind = round(1:(Nt_meas*overlap):(Nt-Nt_meas+1));
t_pred_start_ind = t_meas_start_ind + Nt_meas + Nt_delay;
num_bursts = 1;%length(t_meas_start_ind);

% set least squares model values

%k = logspace(-3,0,51)*2*pi;
k = logspace(-3,0,51)*2*pi;
%k = linspace(1e-3,1e-1,51)*2*pi;
%k = [-k,k];


theta_wavenumber = linspace(-180,180,19)*pi/180;
%theta_wavenumber = linspace(-180,180,19)*pi/180;
%theta_wavenumber = linspace(dir_principal-90,dir_principal+90,9)*pi/180;

theta_wavenumber = theta_wavenumber(2:end);
reg_factor = 5e0;

% for reference
c = sqrt(9.8./k);
%min_time_delay = min(min(circle_radius)./c(2:end));
%max_time_delay = max(max(circle_radius)./c(2:end));

% Set target points
x_target = (round(Nx/4):1:round(Nx/2))*dx;
y_target = (round(Ny/2-20):round(Ny/2+20))*dy;
nx_target = length(x_target);
ny_target = length(y_target);
[x_target,y_target] = meshgrid(x_target,y_target);
x_target = x_target(:);
y_target = y_target(:);
num_points_target = nx_target*ny_target;
x_target_ind = nan(num_points_target,1);
y_target_ind = nan(num_points_target,1);
for i = 1:num_points_target
    x_target_ind(i) = find(Y.x == x_target(i));
    y_target_ind(i) = find(Y.y == y_target(i));
end

% initialize arrays
t_pred = nan(num_bursts,Nt_pred,num_points_target);
z_target_pred = nan(num_bursts,Nt_pred,num_points_target);
z_target_truth = nan(num_bursts,Nt_pred,num_points_target);
x_target_vect = repmat(x_target(:)',[Nt_pred,1]);
y_target_vect = repmat(y_target(:)',[Nt_pred,1]);
time_calc = nan(num_bursts);
% loop over bursts
for i = 1:num_bursts
    % collect buoy array measurements for this burst
    t_target_ind = t_meas_start_ind(i):(t_meas_start_ind(i)+Nt_meas-1);
    x_sentinel_vect = reshape(repmat(x_sentinel(:)',[1,Nt_meas]),[1,Nt_meas*n_points]);
    y_sentinel_vect = reshape(repmat(y_sentinel(:)',[1,Nt_meas]),[1,Nt_meas*n_points]);
    z_sentinel_vect = reshape(z_sentinel(:,t_target_ind),[1,Nt_meas*n_points]);
    t_sentinel_vect = reshape(repmat(Y.t(t_target_ind)',[n_points,1]),[1,Nt_meas*n_points]);
    % set time vector for prediction
    t_pred_ind = t_pred_start_ind(i):(t_pred_start_ind(i)+Nt_pred-1);
    t_pred(i,:,:) = repmat((t_pred_ind'-1)*dt,[1,1,num_points_target]);
    % least squares calculation
    tic
    [z_target_pred_vect,~] = leastSquaresWavePropagation_2D(z_sentinel_vect,...
        t_sentinel_vect,x_sentinel_vect,y_sentinel_vect,...
        t_pred(i,:,:),x_target_vect(:),y_target_vect(:),k,theta_wavenumber,reg_factor);
    z_target_pred(i,:,:) = reshape(z_target_pred_vect,[1,Nt_pred,num_points_target]);
    time_calc(i) = toc;
    % populate ground truth array to compare with prediction
    t_ind_compare = t_pred_ind(t_pred_ind<=Nt & t_pred_ind>=1);
    for j = 1:num_points_target
        z_target_truth(i,t_pred_ind<=Nt & t_pred_ind>=1,j) = squeeze(Y.Z(y_target_ind(j),x_target_ind(j),t_ind_compare));
    end
    %z_circle_check = reshape(z_circle_check,[n_points,Nt]);
end
time_calc_avg = median(time_calc);
%%

target_ind_plot = find(x_target == round(Nx/4)*dx & y_target == round(Ny/2)*dy);
y_target_ind_plot = y_target_ind(target_ind_plot);
x_target_ind_plot = x_target_ind(target_ind_plot);
t_pred_plot = squeeze(t_pred(:,:,target_ind_plot));
z_target_pred_plot = squeeze(z_target_pred(:,:,target_ind_plot));
z_target_truth_plot = squeeze(z_target_truth(:,:,target_ind_plot));
fig(6) = figure(6);
clf(fig(6))
fig(6).Position = [1 1 8 6];
fig(6).PaperPosition = fig(6).Position;
subplot(2,1,1)
hold on
plot(t,squeeze(Y.Z(y_target_ind_plot,x_target_ind_plot,:)))
plot(t_pred_plot,z_target_pred_plot,'-k')
%checkInd = 1;
%plot(t,z_circle(checkInd,:))
%plot(t,z_circle_check(checkInd,:))
hold off
box on
legend('Measured','Predicted')
ylim([-H_sig H_sig])
xlim([min([t';t_pred(:)]),max([t';t_pred(:)])])
xlabel('Time (s)')
ylabel('\eta (m)')

subplot(2,6,7:8)
hold on
pcolor(repmat(Y.x',[Ny,1]),repmat(Y.y,[1,Nx]),Y.Z(:,:,i))
shading flat
%plot(x_target+(-1000:100:1000)*cos(dir_principal*pi/180),...
%    y_target+(-1000:100:1000)*sin(dir_principal*pi/180),'--k')
plot(x_sentinel,y_sentinel,'or')
plot(x_target,y_target,'ok')
hold off
set(gca,'CLim',[-H_sig H_sig],'XLim',[0 Nx*dx],'YLim',[0 Ny*dy])
xlabel('x [m]'), ylabel('y [m]'),

subplot(2,4,7:8)
[z_regress,z_regress_int,z_regress_resid,~,regress_stats] = regress(z_target_truth_plot(:),z_target_pred_plot(:));
regress_r2 = regress_stats(1);
plot(z_target_pred_plot(:),z_target_truth_plot(:),'.')
hold on
plot([-H_sig H_sig],[-H_sig H_sig],'--k')
plot([-H_sig H_sig],[-H_sig H_sig].*z_regress,'--b')
hold off
xlim([-H_sig H_sig])
ylim([-H_sig H_sig])
text(-5*H_sig/6,5*H_sig/6,sprintf('R^2 = %4.1f',100*regress_r2))
xlabel('Predicted')
ylabel('Measured')


print(fig(6),'-dpng',[plotDir '/WEC_Prediction.png'])

%% Compare surface reconstruction
movieLength = Nt_pred; %seconds
for i = 1:movieLength
    fig(1) = figure(1); clf(fig(1));
    fig(1).Position = [2 2 5 10];
    fig(1).PaperPosition = fig(1).Position;
    subplot(2,1,1)
    hold on
    pcolor(reshape(x_target,[ny_target,nx_target]),...
        reshape(y_target,[ny_target,nx_target]),...
        reshape(squeeze(z_target_truth(1,i,:)),[ny_target,nx_target]))
    shading flat
    plot(x_sentinel,y_sentinel,'or')
    hold off
    set(gca,'CLim',1/2*[-H_sig H_sig],'XLim',[0 Nx*dx],'YLim',[0 Ny*dy])
    xlabel('x [m]'), ylabel('y [m]'),
    subplot(2,1,2)
    hold on
    pcolor(reshape(x_target,[ny_target,nx_target]),...
        reshape(y_target,[ny_target,nx_target]),...
        reshape(squeeze(z_target_pred(1,i,:)),[ny_target,nx_target]))
    shading flat
    plot(x_sentinel,y_sentinel,'or')
    hold off
    set(gca,'CLim',1/2*[-H_sig H_sig],'XLim',[0 Nx*dx],'YLim',[0 Ny*dy])
    xlabel('x [m]'), ylabel('y [m]'),
    pause(.01)
end

%% Plot 2D Surface Movie Frames
if plotSurface
    MovieFramesDir = [plotDir '/SurfaceMovieFrames'];
    mkdir(MovieFramesDir)
    movieLength = 30; %seconds
    for i = 1:(movieLength/dt)
        fig(1) = figure(1); clf(fig(1));
        fig(1).Position = [2 2 5 4];
        fig(1).PaperPosition = fig(1).Position;
        hold on
        pcolor(repmat(Y.x',[Ny,1]),repmat(Y.y,[1,Nx]),Y.Z(:,:,i))
        shading flat
        plot(x_target+(-1000:100:1000)*cos(dir_principal*pi/180),...
            y_target+(-1000:100:1000)*sin(dir_principal*pi/180),'--k')
        plot(x_sentinel,y_sentinel,'or')
        plot(x_target,y_target,'ok')
        hold off
        set(gca,'CLim',[-H_sig H_sig],'XLim',[0 Nx*dx],'YLim',[0 Ny*dy])
        xlabel('x [m]'), ylabel('y [m]'),
        print(fig(1),[MovieFramesDir '/SimulationFrame_' sprintf('%03d',i) '.png'],'-dpng')
    end
end

%% Rendered Frames
if plotRenderedFrames
    RenderedMovieFramesDir = [plotDir '/SurfaceMovieFrames_Rendered'];
    mkdir(RenderedMovieFramesDir)
    movieLength = 120; %seconds
    [xs,ys,zs] = sphere;
    
    i1 = 60;
    for i = i1:i1+(movieLength/dt)
        fig(7) = figure(7); clf(fig(7));
        fig(7).Position = [2 2 16 10];
        fig(7).PaperPosition = fig(7).Position;
        subplot(2,1,1)
        hold on
        surf(repmat(Y.x',[Ny,1]),repmat(Y.y,[1,Nx]),Y.Z(:,:,i),Y.Z(:,:,i))
        for j = 1:n_points
            ssphere = surf(x_sentinel(j)+xs,y_sentinel(j)+ys,z_sentinel(j,i)+zs,repmat(shiftdim([1 1 0],-1),[21 21 1]));
        end
        surf(x_target+xs,y_target+ys,Y.Z(y_target_ind,x_target_ind,i)+zs,repmat(shiftdim([1 0 0],-1),[21 21 1]));
        colormap gray
        shading interp
        lighting phong
        material shiny
        %plot(x_circle,y_circle,'or')
        %plot(x_0,y_0,'ok')
        hold off
        %set(gca,'ZLim',[-H_sig H_sig],'CLim',[-H_sig H_sig]);%,'XLim',[200 300])%,'YLim',[200 300])
        axis equal
        set(gca,'ZLim',[-H_sig H_sig],'CLim',[-H_sig H_sig],'XLim',[x_target-300,x_target+50],'YLim',[y_target-150,y_target+200])
        view(0,15)
        %xlabel('x [m]'), ylabel('y [m]'),
        axis off
        %        view(0,5)
        light
        subplot(2,1,2)
        hold on
        plot(t(i1:i),squeeze(Y.Z(y_target_ind,x_target_ind,i1:i)))
        plot(t_pred(t_pred<i),2*z_target_pred(t_pred<i))
        box on
        hold off
        legend('Measured','Predicted')
        xlabel('Time (s)')
        ylabel('\eta [m]')
        ylim([-H_sig H_sig]*2/3)
        %        xlim([min([t';t_pred(:)]),max([t';t_pred(:)])])
        xlim([i1 i1+120])
        print(fig(7),[RenderedMovieFramesDir '/SimulationFrame_' sprintf('%03d',i-i1+1) '.png'],'-dpng')
    end
end