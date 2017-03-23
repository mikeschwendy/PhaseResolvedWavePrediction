clc
clearvars
close all

%% Load data

spectraDir = '/Users/mike/Documents/UW/Research/Data/StereoTests';
allWaverider = load([spectraDir '/Waverider/allDirectionalSpectra.mat']);
load([spectraDir '/SWIFT/allSWIFT2015_withsurfaceelevations.mat']);
stereoDir = '/Volumes/Data/PAPA/TGTcruise2015/MikeScratchTGT/TGT_StereoVideo/StereoXYZResults_FiveMin_WindDirection_BestQC';
resultsDir = '/Users/mike/Documents/UW/Research/Results/SurfaceSimulations';
loadSimulatedData = false;
plotTimeSeries = true;
plotFrequencySpectra = true;
plotStereoComparison = false;
plotLargeDomain = false;
plotRenderedMovie = false;
plotTimeSeriesPrediction = false;

% Enter stereo start time
stereoTime = datenum(2015,1,5,20,30,0); %datenum(2014,12,28,19,29,0); %

% Find data timestamped within 30 minutes of stereo start time
waveriderInds = find(abs(allWaverider.time - stereoTime) < datenum(0,0,0,0,30,0));
swiftInds = find(abs([allSWIFT.time] - stereoTime) < datenum(0,0,0,0,30,0));

% Average waverider directional spectra
EfthetaWaverider = nanmean(allWaverider.Sftheta(:,:,waveriderInds),3);
EfthetaWaverider = pi/180*EfthetaWaverider; % convert to m^2/(Hz*deg)
thetaWaverider = allWaverider.theta;
fWaverider = allWaverider.f;

% Calculate SWIFT directional spectra
if ~isempty(swiftInds)
    [EfthetaSWIFT, thetaSWIFT, fSWIFT, meanDirSWIFT, spreadSWIFT] = ...
        SWIFTdirectionalspectra(allSWIFT(swiftInds));
    thetaSWIFT = wrapTo180(thetaSWIFT);
    [thetaSWIFT, dsort] = sort(thetaSWIFT); % change to -180 to 180
    EfthetaSWIFT = EfthetaSWIFT(:,dsort);
    close all
end
zSwift = allSWIFT(swiftInds(1)).z;
dtSwift = 1/100;

% load stereo data
stereoResultsDir = [stereoDir '/' datestr(stereoTime,'ddmmmyyyy/HHMMUTC') '/XYZResults'];
load([stereoResultsDir '/XYGrid_01.mat'])
zGridFiles = dir([stereoResultsDir '/ZGrid*.mat']);
numMin = numel(zGridFiles);
load([stereoResultsDir '/' char(zGridFiles(1).name)],'zgrid')
[NyStereo,NxStereo,NtStereo] = size(zgrid);
fpsStereo = NtStereo/60;
zgridTot = nan(NyStereo,NxStereo,fpsStereo*60*numMin);
zgridTot(:,:,1:NtStereo) = zgrid;
for j = 2:numMin
    load([stereoResultsDir '/' char(zGridFiles(j).name)],'zgrid')
    zgridTot(:,:,((j-1)*NtStereo+1):(j*NtStereo)) = zgrid;
end
dtStereo = 1/fpsStereo;
tStereo = (1:(numMin*NtStereo))'*dtStereo;

mkdir([resultsDir '/' datestr(stereoTime,'ddmmmyyyy')])
plotDir = [resultsDir '/' datestr(stereoTime,'ddmmmyyyy/HHMMUTC')];
mkdir(plotDir)

%% Transform into WAFO spec struct format
S_Waverider = struct();
S_Waverider.date = stereoTime;
S_Waverider.type = 'dir';
S_Waverider.S = EfthetaWaverider'*180/pi;
S_Waverider.f = fWaverider;
S_Waverider.theta = thetaWaverider'*pi/180;
S_Waverider.phi = 0;

S_Swift = struct();
S_Swift.date = stereoTime;
S_Swift.type = 'dir';
S_Swift.S = EfthetaSWIFT'*180/pi;
S_Swift.f = fSWIFT;
S_Swift.theta = thetaSWIFT'*pi/180;
S_Swift.phi = 0;

meanDir = dspec2char(S_Waverider,'Mdir');
meanDir = cell2mat(meanDir);
% theta = 0 --> wave propogates in -x,
% theta = pi/2 --> wave propogates in -y


S=jonswap;
D=spreading(linspace(-pi,pi,51),'cos2s');
S=mkdspec(S,D,0);
%Srot = rotspec(S,meanDir,1);
Srot = rotspec(S,pi/2,1);
S_Waverider = ttspec(S);
Srot = ttspec(Srot);

% Srot.S(S.theta<0.95*pi/2 | S.theta>1.05*pi/2,:) = 0;

%% Plot Directional Spectra
fig(4) = figure(4); clf(fig(4));
fig(4).Position = [2 2 16 8];
fig(4).PaperPosition = fig(4).Position;
s1 = subplot(1,2,1);
polarPcolor(S_Waverider.f',S_Waverider.theta'*180/pi,(S_Waverider.S));
s2 = subplot(1,2,2);
polarPcolor(Srot.f',Srot.theta'*180/pi,(Srot.S));
print(fig(4),'-dpng',[plotDir '/RotatedDirectionalSpectrum.png'],'-r300')

S_Waverider = Srot;

%% Define simulation domain
wMax = max(2*pi*S_Waverider.f);
wMin = min(2*pi*S_Waverider.f);
kMax = w2k(wMax);
kMin = w2k(wMin);
lMin = 2*pi/kMax;  % dx should equal approx. lMin/4
dx = 1;
lMax = 2*pi/kMin;
L = 200;
Nx = round(L/dx);
Ny = Nx;
dy = dx;
tMin = 2*pi/wMax;  % dt should equal approx. tMin/4
dt = 1;
fpsSim = 1/dt;
tMax = 2*pi/wMin;
T = 1800;
Nt = round(T/dt);

%% Call WAFO codes to simulate 3-D surface

if ~loadSimulatedData
fftdim = 1; % if fftdim = 2, 2D ifft at each time.  if fftdim = 1, 1D ifft at each space
[Y,~] = seasim(S_Waverider,Nx,Ny,Nt,dx,dy,dt,fftdim,0);
save([plotDir '/SimulationData.mat'],'Y')
else
    load([plotDir '/SimulationData.mat'],'Y')
end
%% Check Sample Time Series
if plotTimeSeries
    fig(1) = figure(1); clf(fig(1));
    s1(1) = subplot(3,1,1);
    plot(Y.t,squeeze(Y.Z(ceil(Nx/2),ceil(Ny/2),:)),'-k')
    s1(1).XLim = [0 numMin*60];
    s1(1).YLim = [-4 4];
    ylabel('z [m]')
    title('Simulation')
    s1(2) = subplot(3,1,2);
    plot(tStereo,squeeze(zgridTot(ceil(NyStereo/2),ceil(NxStereo/2),:)),'-b')
    s1(2).XLim = s1(1).XLim;
    s1(2).YLim = s1(1).YLim;
    ylabel('z [m]')
    title('Stereo')
    s1(3) = subplot(3,1,3);
    plot((1:length(zSwift))*dtSwift,zSwift,'-r')
    s1(3).XLim = s1(1).XLim;
    s1(3).YLim = s1(1).YLim;
    ylabel('z [m]')
    title('SWIFT')
    xlabel('t [s]')
    print(fig(1),[plotDir '/SampleTimeSeries.png'],'-dpng','-r300')
end
%% Check Frequency Spectra
binAvg = 1;
windowLengthSeconds = 120;
windowLengthSim = 2^nextpow2(windowLengthSeconds*fpsSim);
fSim = 0:(binAvg*fpsSim/windowLengthSim):(fpsSim/2);
nfSim = length(fSim);
xSubSim = round(linspace(1,Nx,min(40,Nx)));
ySubSim = round(linspace(1,Ny,min(40,Ny)));
SfSim = nan(length(ySubSim),length(xSubSim),nfSim);
for j = 1:length(xSubSim)
    for k = 1:length(ySubSim)
        %SfSim2(k,j) = dat2spec([Y.t,squeeze(Y.Z(ySubSim(k),xSubSim(j),:))],128,[],[],[],'psdo',[],'f');
        % 2*pwelch because giving frequency vector requires two-sided
        % psd
        SfSim(k,j,:) = 2*pwelch(squeeze(Y.Z(ySubSim(k),xSubSim(j),:)),windowLengthSim,0.5*windowLengthSim,fSim,fpsSim,'psd');
    end
end

windowLengthStereo = 2^nextpow2(windowLengthSeconds*fpsStereo);
fStereo = 0:(binAvg*fpsStereo/windowLengthStereo):(fpsStereo/2);
nfStereo = length(fStereo);
xSubStereo = round(linspace(1,NxStereo,min(40,NxStereo)));%1:4:NxStereo;
ySubStereo = round(linspace(1,NyStereo,min(40,NyStereo)));%1:4:NyStereo;
SfStereo = nan(length(ySubStereo),length(xSubStereo),nfStereo);
for j = 1:length(xSubStereo)
    for k = 1:length(ySubStereo)
        SfStereo(k,j,:) = 2*pwelch(squeeze(zgridTot(ySubStereo(k),xSubStereo(j),:)),windowLengthStereo,0.5*windowLengthStereo,fStereo,fpsStereo,'psd');
    end
end
Sf_Waverider = spec2spec(S_Waverider,'freq');
Sf_Swift = spec2spec(S_Swift,'freq');
HsWaverider0 = 4*sqrt(trapz(Sf_Waverider.f,Sf_Waverider.S))
HsSim1 = 4*nanstd(Y.Z(:))
HsSim2 = 4*sqrt(trapz(fSim,squeeze(nanmean(nanmean(SfSim,2),1))))
HsStereo1 = 4*nanstd(zgridTot(:))
HsStereo2 = 4*sqrt(trapz(fStereo,squeeze(nanmean(nanmean(SfStereo,2),1))))
HsSwift1 = 4*nanstd(zSwift)
HsSwift2 = 4*sqrt(trapz(Sf_Swift.f,Sf_Swift.S))
%% Plot Frequency Spectra
if plotFrequencySpectra
    fig(2) = figure(2); clf(fig(2));
    plot(Sf_Waverider.f,Sf_Waverider.S,'-k','linewidth',2)
    hold on
    plot(fSim,squeeze(nanmean(nanmean(SfSim,2),1)),'-k')
    %plot(fStereo,squeeze(nanmean(nanmean(SfStereo,2),1)),'-b')
    %plot(Sf_Swift.f,Sf_Swift.S,'-r')
    hold off
    %set(gca,'YScale','log','XScale','log','xlim',[0.05,0.6],'ylim',[1e-4 1e3])
    set(gca,'YScale','linear','XScale','linear','xlim',[0,0.3]);%,'ylim',[0 30])
    legend('Original Waverider Spectrum','Simulation','Stereo','SWIFT',...
        'location','best')
    print(fig(2),[plotDir '/CheckFrequencySpectra.png'],'-dpng','-r300')
end
%% Surface Movie Frames - Comparison with Stereo

if plotStereoComparison
    StereoMovieFramesDir = [plotDir '/2DSurfaceMovieFrames_StereoComparison'];
    mkdir(StereoMovieFramesDir)
    %tPlot = dt:dt:20;
    for j = 1:fpsSim*30
        jSim = j;%find(Y.t == tPlot(j));
        jStereo = j;%find(tStereo == tPlot(j));
        fig(3) = figure(3); clf(fig(3));
        fig(3).Position = [2 2 8 4];
        fig(3).PaperPosition = fig(3).Position;
        subplot(1,2,1)
        pcolor(xgrid,ygrid,zgridTot(:,:,j))
        shading('flat')
        set(gca,'CLim',[-4 4],'XLim',[-15 15],'YLim',[-20 20])
        xlabel('x [m]'), ylabel('y [m]'), title('Stereo Data')
        subplot(1,2,2)
        pcolor(Y.x,Y.y,Y.Z(:,:,j))
        shading('flat')
        set(gca,'CLim',[-4 4],'XLim',[5 35],'YLim',[0 40])
        xlabel('x [m]'), ylabel('y [m]'), title('Simulation')
        print(fig(3),[StereoMovieFramesDir '/SimulationStereoComparison_' sprintf('%02d',j) '.png'],'-dpng')
    end
end

%% Surface Movie Frames - Large Domain
if plotLargeDomain
    LargeMovieFramesDir = [plotDir '/2DSurfaceMovieFrames_LargeDomain'];
    mkdir(LargeMovieFramesDir)
    for j = 1:fpsSim*30
        jSim = j;%find(Y.t == tPlot(j));
        fig(5) = figure(5); clf(fig(5));
        fig(5).Position = [2 2 5 4];
        fig(5).PaperPosition = fig(5).Position;
        pcolor(Y.x,Y.y,Y.Z(:,:,j))
        shading('flat')
        set(gca,'CLim',[-4 4],'XLim',[0 Nx*dx],'YLim',[0 Ny*dy])
        cbar = colorbar('location','northoutside');
        xlabel('x [m]'), ylabel('y [m]'), %title('Simulation')
        xlabel(cbar,'z [m]')
        axis equal
        print(fig(5),[LargeMovieFramesDir '/SimulationLargeDomain_' sprintf('%02d',j) '.png'],'-dpng')
    end
end

%%
if plotRenderedMovie
    RenderedFramesDir = [plotDir '/2DSurfaceMovieFrames_Rendered'];
    mkdir(RenderedFramesDir)
    for j = 1:fpsSim*30
        fig(7) = figure(7); clf(fig(7));
        fig(7).Position = [2 2 5 4];
        fig(7).PaperPosition = fig(7).Position;
        surf(Y.x,Y.y,Y.Z(:,:,j))
        hold on
        [xs,ys,zs] = sphere;
        xInd = find(Y.x == 10);
        yInd = find(Y.y == 100);
        ssphere = surf(5*xs+Y.x(xInd),5*ys+Y.y(yInd),5*zs+Y.Z(yInd,xInd,j),repmat(shiftdim([1 1 0],-1),[21 21 1]));
        colormap gray
        shading interp
        lighting phong
        material shiny
        view(0,25)
        set(gca,'clim',[-5 5],'xgrid','off','ygrid','off','zgrid','off')
        axis equal
        axis off
        axis([0 L 0 L -10 10])
        light
        hold off
        %pause
        print(fig(7),[RenderedFramesDir '/SimulationRendered_' sprintf('%02d',j) '.png'],'-dpng')
    end
end

%%
if plotTimeSeriesPrediction
    fig(8) = figure(8); clf(fig(8));
    fig(8).Position = [2 2 12 4];
    fig(8).PaperPosition = fig(8).Position;
    xInd = find(Y.x == 10);
    xInd2 = find(Y.x == 60);
    yInd = find(Y.y == 100);
    yInd2 = find(Y.y == 100);
    yInd3 = find(Y.y == 120);
    yInd4 = find(Y.y == 80);
    plot(Y.t,squeeze(Y.Z(yInd,xInd,:)),'-k','linewidth',2);
    %         hold on
    %         plot(Y.t,squeeze(Y.Z(yInd2,xInd2,:)),'-b');
    %         plot(Y.t,squeeze(Y.Z(yInd3,xInd2,:)),'-r');
    %         plot(Y.t,squeeze(Y.Z(yInd4,xInd2,:)),'-g');
    %         hold off
    set(gca,'ylim',[-3 3],'xlim',[0 60])
    xlabel('Time (s)')
    ylabel('Z (m)')
    hold off
    print(fig(8),[plotDir '/PredictionTimeSeries.png'],'-dpng','-r300')
end

%% Calculate directional spectrum
W = [Y.t,reshape(Y.Z,[Nx*Ny,Nt])'];
pos = [reshape(repmat(Y.x',[Ny,1]),[Nx*Ny,1]),reshape(repmat(Y.y,[1,Nx]),[Nx*Ny,1]),zeros(Nx*Ny,1),ones(Nx*Ny,1),ones(Nx*Ny,1)];
h = inf;
windowLengthSeconds = 120;
windowLengthSim = 2^nextpow2(windowLengthSeconds*fpsSim);
Nfft = 2*windowLengthSim;
Ntheta = 90;
method = 'MLM';
opt = struct('window',2*windowLengthSim);
[S_test,D_test,Sw_test,Fcof_test] = dat2dspec(W,pos,h,Nfft,Ntheta,method);%,opt);
S_test = ttspec(S_test,'f');
%%
fig(6) = figure(6); clf(fig(6));
fig(6).Position = [2 2 16 8];
fig(6).PaperPosition = fig(6).Position;
s1 = subplot(1,2,1);
polarPcolor(S_Waverider.f(S_Waverider.f<0.3)',S_Waverider.theta'*180/pi,log10(S_Waverider.S(:,S_Waverider.f<0.3)));
set(gca,'clim',[0 2])
s2 = subplot(1,2,2);
polarPcolor(S_test.f(S_test.f<0.3)',S_test.theta'*180/pi,log10(S_test.S(:,S_test.f<0.3)));
set(gca,'clim',[0 2])
print(fig(6),'-dpng',[plotDir '/RecalculatedDirectionalSpectrum.png'],'-r300')

%%
 S  = jonswap;
 D  = spreading(linspace(-pi,pi,51),'cos2s');
 Sd = mkdspec(S,D,0);
 Nx = 3; Ny = 2; Nt = 2^15; dx = 10; dy = 10;dt = 0.5;
 F  = seasim(Sd,Nx,Ny,Nt,dx,dy,dt,1,0);  
 Z  = permute(F.Z,[3 1 2]);
 [X,Y] = meshgrid(F.x,F.y);
 N = Nx*Ny;
 types = repmat(sensortypeid('n'),N,1);
 bfs   = ones(N,1);
 pos   = [X(:),Y(:),zeros(N,1)];
 h = inf;
 Se = dat2dspec([F.t Z(:,:)],[pos types,bfs],h,256,101); % seasim is possibly wrong

fig(6) = figure(6); clf(fig(6));
fig(6).Position = [2 2 16 8];
fig(6).PaperPosition = fig(6).Position;
wMax = pi/2;
s1 = subplot(1,2,1);
polarPcolor(Sd.w(Sd.w<wMax)',Sd.theta'*180/pi,log10(Sd.S(:,Sd.w<wMax)));
set(gca,'clim',[0 2])
s2 = subplot(1,2,2);
polarPcolor(Se.w(Se.w<wMax)',Se.theta'*180/pi,log10(Se.S(:,Se.w<wMax)));
set(gca,'clim',[0 2])
