clc
clearvars
close all


%% choose data source as SWIFT, CDIP, or Waverider
% dataSource = 'CDIP';
% sourceDirectory = '/Users/mike/Dropbox/SanNicolasIsland/CDIP 138 Files';
% resultsDirectory = '/Users/mike/Documents/UW/Research/Results/SurfaceSimulations_2D/CDIP_SanNicolasIsland_Mar2016';

% dataSource = 'SWIFT';
% sourceDirectory = '/Users/mike/Documents/UW/Research/Results/PapaSwiftVisualization';
% resultsDirectory = '/Users/mike/Documents/UW/Research/Results/SurfaceSimulations_2D/SWIFT_StationPapa_Jan2015';

dataSource = 'Waverider';
sourceDirectory = '/Users/mike/Documents/UW/Research/Data/StereoTests/Waverider';
resultsDirectory = '/Users/mike/Documents/UW/Research/Results/SurfaceSimulations_2D/Waverider_StationPapa_Jan2015';

%% Load spectra and put in wafo format
switch dataSource
    case 'Waverider'
        % Load Waverider
        allWaverider = load([sourceDirectory '/allDirectionalSpectra.mat']);
        
        % Enter start time
        startTime = datenum(2015,1,5,20,30,0); %datenum(2014,12,28,19,29,0); %
        
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
        S_Measured.f = fWaverider;
        S_Measured.theta = thetaWaverider'*pi/180;
        S_Measured.phi = 0;
        
    case 'SWIFT'
        % Load SWIFT
        swiftInd = '01052';
        load([sourceDirectory '/SWIFTdata_' swiftInd '.mat'],'SWIFT1');
        
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
        cdipFile = [sourceDirectory '/March_2016.dat'];
        % date/time string MUST be of the form 'YYYYMMDDHHmm'
        % cdipTimestring = '201603010121';
        % cdipTimestring = '201603020251';
        % cdipTimestring = '201603030121';
        % cdipTimestring = '201603040121';
        % cdipTimestring = '201603040151';
        cdipTimestring = '201603041521';
        [tspc,tsys] = read_cdip_buoy(cdipFile,cdipTimestring);
        [~, Etheta] = MEM_directionalestimator(tspc.a1,tspc.a2,tspc.b1,tspc.b2,tspc.en,0);
        dtheta = 2;
        theta = -(-180:dtheta:179);  % start with cartesion (a1 is positive east velocities, b1 is positive north)
        
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

%% Spectral parameters

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

% For reference...
wMax = max(2*pi*S_Measured.f);
wMin = min(2*pi*S_Measured.f);
kMax = w2k(wMax);
kMin = w2k(wMin);
% dx need not be below lMin/4
lMin = 2*pi/kMax;
lMax = 2*pi/kMin;
% dt need not be below tMin/4
tMin = 2*pi/wMax;
tMax = 2*pi/wMin;

%% Perform Simulation
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

%% Save simulation results
% waverider output directory
outputDirectory = [resultsDirectory '/' dataSource '_' datestr(S_Measured.date,'ddmmmyyyy') '_' sprintf('L%d_dx%d_T%d_dt%d',L,dx,T,dt)];
mkdir(outputDirectory)
save([outputDirectory '/SimulationData.mat'],'Y','dx','dt','L','T','Nx','Nt','x','t','Ny','dy','y','S_Measured')

%% Plot directional spectra
fig(1) = figure(1); clf(fig(1));
pcolor(S_Measured.theta*180/pi,S_Measured.f,log10(S_Measured.S'))
%pcolor(S_Measured.theta*180/pi,S_Measured.f,S_Measured.S')
shading('flat')
hold on
plot(dir_peak,f_peak,'xk')
text(dir_peak+5,f_peak,'Peak','color','k')
plot(dir_mean,f_mean,'xr')
text(dir_mean+5,f_mean,'Mean','color','r')
hold off
xlabel('\theta (degrees)')
ylabel('f (Hz)')
print(fig(1),[outputDirectory '/DirectionalSpectra.png'],'-dpng')

%% Check Frequency Spectra
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
%    plot(Sf_Measured.f,nansum(S_Measured.S*dtheta),'-k','linewidth',2)
plot(Sf_Measured.f,Sf_Measured.S,'-k','linewidth',2)
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
print(fig(2),[outputDirectory '/FrequencySpectra.png'],'-dpng','-r300')
