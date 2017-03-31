clc
clearvars
close all

plotSurface = false;
plotRenderedFrames = false;

% Load simulation data
SimulationDirectory = '/Users/mike/Documents/UW/Research/Results/SurfaceSimulations_2D/Waverider_StationPapa_Jan2015/Waverider_05Jan2015_L500_dx4_T300_dt1';
load([SimulationDirectory '/SimulationData.mat'],'Y','dx','dt','L','T','Nx','Nt','x','t','Ny','dy','y','S_Measured')

% Choose hypothetical buoy locations
arrayShape = 'circle';
%arrayShape = 'box';

switch arrayShape
    case 'circle'
    % Setup buoy array - Circle
    x_center_ind = round(Nx/2);
    y_center_ind = round(Ny/2);
    x_center = x(x_center_ind);
    y_center = y(y_center_ind);
    n_theta = 5;
    theta_spread = 120; % deg
    dir_principal = cell2mat(dspec2char(S_Measured,'Wdir'))*180/pi; % peak direction
    theta = linspace(dir_principal-theta_spread/2,dir_principal+theta_spread/2,n_theta)';
    circle_radius = [100, 150]; % m
    [x_array,y_array,z_array] = setBuoyArray_Circle(x_center,y_center,theta,circle_radius,dx,dy,Y);

    case 'box'
    % Setup buoy array - box
    dir_principal = cell2mat(dspec2char(S_Measured,'Wdir'))*180/pi; % peak direction
    box_side_length = 200;
    x_center_ind = round(Nx/2);
    y_center_ind = round(Ny/2);
    x_center = x(x_center_ind);
    y_center = y(y_center_ind);
    [x_array,y_array,z_array] = setBuoyArray_Box(x_center,y_center,box_side_length,num_buoys,dir_principal,dx,dy,Y);

end


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
n_points = length(x_array);
% loop over bursts
for i = 1:num_bursts
    % collect buoy array measurements for this burst
    t_target_ind = t_meas_start_ind(i):(t_meas_start_ind(i)+Nt_meas-1);
    x_array_vect = reshape(repmat(x_array(:)',[1,Nt_meas]),[1,Nt_meas*n_points]);
    y_array_vect = reshape(repmat(y_array(:)',[1,Nt_meas]),[1,Nt_meas*n_points]);
    z_array_vect = reshape(z_array(:,t_target_ind),[1,Nt_meas*n_points]);
    t_array_vect = reshape(repmat(Y.t(t_target_ind)',[n_points,1]),[1,Nt_meas*n_points]);
    % set time vector for prediction
    t_pred_ind = t_pred_start_ind(i):(t_pred_start_ind(i)+Nt_pred-1);
    t_pred(i,:,:) = repmat((t_pred_ind'-1)*dt,[1,1,num_points_target]);
    % least squares calculation
    tic
    [z_target_pred_vect,~] = leastSquaresWavePropagation_2D(z_array_vect,...
        t_array_vect,x_array_vect,y_array_vect,...
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
plot(x_array,y_array,'or')
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
    plot(x_array,y_array,'or')
    hold off
    set(gca,'CLim',1/2*[-H_sig H_sig],'XLim',[0 Nx*dx],'YLim',[0 Ny*dy])
    xlabel('x [m]'), ylabel('y [m]'),
    subplot(2,1,2)
    hold on
    pcolor(reshape(x_target,[ny_target,nx_target]),...
        reshape(y_target,[ny_target,nx_target]),...
        reshape(squeeze(z_target_pred(1,i,:)),[ny_target,nx_target]))
    shading flat
    plot(x_array,y_array,'or')
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
        plot(x_array,y_array,'or')
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
            ssphere = surf(x_array(j)+xs,y_array(j)+ys,z_array(j,i)+zs,repmat(shiftdim([1 1 0],-1),[21 21 1]));
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