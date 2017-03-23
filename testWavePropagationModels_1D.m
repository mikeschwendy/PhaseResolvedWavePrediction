clc
clearvars
close all

% Load waverider data from Station Papa Jan. 2015
spectraDir = '/Users/mike/Documents/UW/Research/Data/StereoTests';
allWaverider = load([spectraDir '/Waverider/allDirectionalSpectra.mat']);
loadSimulatedData = true;
plotSurface = false;
plotFrequencySpectra = true;
plotFFTMovie = false;
resultsDir = '/Users/mike/Documents/UW/Research/Results/SurfaceSimulations_1D';

% Enter start time
startTime = datenum(2015,1,5,20,30,0); %datenum(2014,12,28,19,29,0); %
mkdir([resultsDir '/' datestr(startTime,'ddmmmyyyy')])
plotDir = [resultsDir '/' datestr(startTime,'ddmmmyyyy/HHMMUTC')];
mkdir(plotDir)

% Find data timestamped within 30 minutes of start time
waveriderInds = find(abs(allWaverider.time - startTime) < datenum(0,0,0,0,30,0));

% Average waverider directional spectra
EfthetaWaverider = nanmean(allWaverider.Sftheta(:,:,waveriderInds),3);
EfthetaWaverider = pi/180*EfthetaWaverider; % convert to m^2/(Hz*deg)
thetaWaverider = allWaverider.theta;
fWaverider = allWaverider.f;

%% Transform into WAFO spec struct format
S_Waverider = struct();
S_Waverider.date = startTime;
S_Waverider.type = 'dir';
%S_Waverider.S = EfthetaWaverider'*180/pi;
S_Waverider.S = zeros(length(thetaWaverider),length(fWaverider));
S_Waverider.S((length(thetaWaverider)+1)/2,:) = sum(pi/180*4*EfthetaWaverider,2)'*180/pi;
%S_Waverider.S(1,:) = sum(pi/180*4*EfthetaWaverider,2)'*180/pi;
S_Waverider.f = fWaverider;
S_Waverider.theta = thetaWaverider'*pi/180;
S_Waverider.phi = 0;

%% Define simulation domain
wMax = max(2*pi*S_Waverider.f);
wMin = min(2*pi*S_Waverider.f);
kMax = w2k(wMax);
kMin = w2k(wMin);
lMin = 2*pi/kMax;  % dx should equal approx. lMin/4
lMax = 2*pi/kMin;
tMin = 2*pi/wMax;  % dt should equal approx. tMin/4
tMax = 2*pi/wMin;
dx = 1;
dt = 0.2;
fpsSim = 1/dt;
L = 500;
T = 1800;
Nx = round(L/dx);
Nt = round(T/dt);
x = linspace(0,L-dx,Nx); % m
t = linspace(0,T-dt,Nt); % s

%% Call WAFO codes to simulate surface trace
if ~loadSimulatedData
    fftdim = 1; % if fftdim = 2, 2D ifft at each time.  if fftdim = 1, 1D ifft at each space
    [Y,~] = seasim(S_Waverider,Nx,0,Nt,dx,0,dt,fftdim,0);
    save([plotDir '/SimulationData.mat'],'Y')
else
    load([plotDir '/SimulationData.mat'],'Y')
end

%% Check Frequency Spectra
binAvg = 1;
windowLengthSeconds = 120;
windowLengthSim = 2^nextpow2(windowLengthSeconds*fpsSim);
fSim = 0:(binAvg*fpsSim/windowLengthSim):(fpsSim/2);
nfSim = length(fSim);
xSubSim = round(linspace(1,Nx,min(40,Nx)));
SfSim = nan(length(xSubSim),nfSim);
for j = 1:length(xSubSim)
        % 2*pwelch because giving frequency vector requires two-sided
        SfSim(j,:) = 2*pwelch(squeeze(Y.Z(xSubSim(j),:)),windowLengthSim,0.5*windowLengthSim,fSim,fpsSim,'psd');
end
Sf_Waverider = spec2spec(S_Waverider,'freq');

%% Plot Frequency Spectra
if plotFrequencySpectra
    fig(2) = figure(2); clf(fig(2));
    %plot(Sf_Waverider.f,Sf_Waverider.S,'-k','linewidth',2)
    hold on
    plot(Sf_Waverider.f,nansum(S_Waverider.S*4*pi/180),'-k','linewidth',2)
    plot(fSim,squeeze(nanmean(SfSim,1)),'-k')
    hold off
    %set(gca,'YScale','log','XScale','log','xlim',[0.05,0.6],'ylim',[1e-4 1e3])
    set(gca,'YScale','linear','XScale','linear','xlim',[0,0.3]);%,'ylim',[0 30])
    legend('Original Waverider Spectrum','Simulation','Stereo','SWIFT',...
        'location','best')
    %print(fig(2),[plotDir '/CheckFrequencySpectra.png'],'-dpng','-r300')
end
fMean = sum(squeeze(nanmean(SfSim,1)).*fSim)/sum(squeeze(nanmean(SfSim,1)));
kMean = (2*pi*fMean).^2/9.8;
cMean = 2*pi*fMean/kMean;
%% Plot 1D Surface Movie Frames
if plotSurface
    MovieFramesDir = [plotDir '/SurfaceMovieFrames'];
    mkdir(MovieFramesDir)
    movieLength = 30; %seconds
    for i = 1:(fpsSim*movieLength)
        fig(1) = figure(1); clf(fig(1));
        fig(1).Position = [2 2 10 4];
        fig(1).PaperPosition = fig(1).Position;
        hold on
        plot(Y.x,Y.Z(:,i))
        plot(0+i*cMean/fpsSim,3,'x')
        hold off
        set(gca,'YLim',[-4 4],'XLim',[0 Nx*dx])
        xlabel('x [m]'), ylabel('z [m]'),
        %axis equal
        print(fig(1),[MovieFramesDir '/SimulationFrame_' sprintf('%03d',i) '.png'],'-dpng')
    end
end

%% Test linear propagation model (i.e. Brown and Jensen 2001)

f = fpsSim*([0:(Nt/2-1),-Nt/2:-1])/Nt;  % 1/s
k = sign(f).*(2*pi*f).^2/(9.8);
x0 = Nx/2;
A = fft(Y.Z(x0,:));
Z_prop = ifft(repmat(A,[Nx,1]).*exp(1i*(repmat(x',[1,Nt])-x(x0)).*repmat(k,[Nx,1])),[],2);

%%
if plotFFTMovie 
movieLength = T; %seconds
for i = 1:(fpsSim*movieLength)
    fig(1) = figure(1); clf(fig(1));
    fig(1).Position = [2 2 10 4];
    fig(1).PaperPosition = fig(1).Position;
    hold on
    plot(Y.x,Y.Z(:,i))
    plot(Y.x,real(Z_prop(:,i)),'--')
    plot(mod(-i*cMean/fpsSim,Nx*dx),3,'x')
    hold off
    set(gca,'YLim',[-4 4],'XLim',[0 Nx*dx])
    xlabel('x [m]'), ylabel('z [m]'),
    pause(.01)
    %axis equal
    %print(fig(1),[MovieFramesDir '/SimulationFrame_' sprintf('%03d',i) '.png'],'-dpng')
end
end
%% for xSep = 0:50
xSep = 150;
fig(3) = figure(3); 
fig(3).Position = [2 2 14 10];
fig(3).PaperPosition = fig(3).Position;
clf(fig(3))
subplot(2,2,1)
hold on
plot(t,Y.Z(x0+xSep,:))
plot(t,real(Z_prop(x0+xSep,:)),'r--')
hold off
%xlim([dt*Nt/2 dt*Nt/2+60])
ylim([-4 4])
xlabel('t [s]')
ylabel('\eta [m]')
subplot(2,2,3)
hold on
plot(Y.Z(x0+xSep,:),real(Z_prop(x0+xSep,:)),'.')
plot([-4 4],[-4 4],'--k')
hold off
set(gca,'xlim',[-4 4],'ylim',[-4 4])
xlabel('\eta (actual)')
ylabel('\eta (predicted)')

tCheck = 8500;
subplot(2,2,2)
hold on
plot(x,Y.Z(:,tCheck))
plot(x,real(Z_prop(:,tCheck)),'r--')
hold off
%xlim([dt*Nt/2 dt*Nt/2+60])
ylim([-4 4])
xlabel('x [m]')
ylabel('\eta [m]')
subplot(2,2,4)
hold on
plot(Y.Z(:,tCheck),real(Z_prop(:,tCheck)),'.')
plot([-4 4],[-4 4],'--k')
hold off
set(gca,'xlim',[-4 4],'ylim',[-4 4])
xlabel('\eta (actual)')
ylabel('\eta (predicted)')

%% Least squares calculation
clc
%f_ls = [fSim(1:100)+2e-3,fSim(1:100)];
%k_ls = [-(2*pi*f_ls(1:100)).^2/(9.8),(2*pi*f_ls(101:200)).^2/(9.8)];

f_ls = [-fSim(11:2:100)+1e-3,fSim(10:2:100)];
k_ls = -sign(f_ls).*(2*pi*f_ls).^2/(9.8);


%f_ls = [fSim(1:100)];
%k_ls = -sign(f_ls).*(2*pi*f_ls).^2/(9.8);

%f_ls = [fSim(2:100)];
%k_ls = [-sign(f_ls).*(2*pi*f_ls).^2/(9.8),+sign(f_ls).*(2*pi*f_ls).^2/(9.8)];
%f_ls = [f_ls,f_ls];

%k_ls = -25:0.125:25;
%f_ls = sqrt(9.8*abs(k_ls))/(2*pi);


%f_ls = -fSim(1:200);
%k_ls = [-(2*pi*f_ls).^2/(9.8), (2*pi*(f_ls+1e-3)).^2/(9.8)];
%f_ls = [f_ls,f_ls+1e-3];
Nf_ls = length(f_ls);
[t_tot,x_tot] = meshgrid(t,x);

x_meas_ind = [Nx/2,Nx/2+10];%(Nx-100):10:Nx;
Nx_meas = length(x_meas_ind);
t_meas_ind = (Nt-2000):(Nt-1500);
Nt_meas = length(t_meas_ind);
z_meas = Y.Z(x_meas_ind,t_meas_ind);
x_meas = x_tot(x_meas_ind,t_meas_ind);
t_meas = t_tot(x_meas_ind,t_meas_ind);

k_meas_tot = repmat(k_ls,[Nx_meas*Nt_meas,1]);
f_meas_tot = repmat(f_ls,[Nx_meas*Nt_meas,1]);
t_meas_tot = repmat(t_meas(:),[1,Nf_ls]);
x_meas_tot = repmat(x_meas(:),[1,Nf_ls]);
P_meas = exp(1i.*(x_meas_tot.*k_meas_tot-t_meas_tot.*f_meas_tot*2*pi));
A = P_meas\z_meas(:); 
%A = pinv(conj(P_meas')*P_meas)*conj(P_meas')*z_meas(:);
%A = pinv(A)*z_meas(:);
z_check = reshape(real(P_meas*A),[Nx_meas,Nt_meas]);

xSep2 = -10;
tSep2 = 0;%round((xSep2*dx/cMean)/dt);
x_pred_ind = x_meas_ind(1)+xSep2;%[x_meas_ind-xSep2,x_meas_ind+xSep2];
Nx_pred = length(x_pred_ind);
t_pred_ind = t_meas_ind+tSep2;
Nt_pred = length(t_pred_ind);
z_actual = Y.Z(x_pred_ind,t_pred_ind);
x_pred = x_tot(x_pred_ind,t_pred_ind);
t_pred = t_tot(x_pred_ind,t_pred_ind);


k_pred_tot = repmat(k_ls,[Nx_pred*Nt_pred,1]);
f_pred_tot = repmat(f_ls,[Nx_pred*Nt_pred,1]);
t_pred_tot = repmat(t_pred(:),[1,Nf_ls]);
x_pred_tot = repmat(x_pred(:),[1,Nf_ls]);
P_pred = exp(1i.*(x_pred_tot.*k_pred_tot-t_pred_tot.*f_pred_tot*2*pi));
z_pred = reshape(real(P_pred*A),[Nx_pred,Nt_pred]);

fig(6) = figure(6); 
clf(fig(6))
subplot(2,1,1)
hold on
plot(t_meas(1,:)',z_meas(1,:)')
plot(t_meas(1,:)',z_check(1,:)')
plot(t_pred',z_actual')
plot(t_pred',z_pred')

% plot(x_meas(:,1)',z_meas(:,1)')
% plot(x_meas(:,1)',z_check(:,1)')
% plot(x_pred',z_actual')
% plot(x_pred',z_pred')
hold off
legend('Measured','Check','Target, actual','Target, prediction')
ylim([-2 2])
xlim([1400 1460])
subplot(2,1,2)
plot(t_pred',abs(z_actual'-z_pred'))
%plot(z_pred'./z_actual','.')
%plot(z_meas(1,:)'./z_check(1,:)','.')
ylim([0 1])