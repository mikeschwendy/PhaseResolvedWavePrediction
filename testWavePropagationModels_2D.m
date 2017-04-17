clc
clearvars
close all

plotPointComparison = true;
plotSurfaceComparison = false;
plotRenderedFrames = false;

% Load simulation data
%SimulationDirectory = '/Users/mike/Documents/UW/Research/Results/SurfaceSimulations_2D/Waverider_StationPapa_Jan2015/Waverider_05Jan2015_L500_dx4_T300_dt1';
SimulationDirectory = '/Users/mike/Documents/UW/Research/Results/SurfaceSimulations_2D/CDIP_SanNicolasIsland_Mar2016/CDIP_04Mar2016_L500_dx4_T300_dt1';
load([SimulationDirectory '/SimulationData.mat'],'Y','dx','dt','L','T','Nx','Nt','x','t','Ny','dy','y','S_Measured')
plotDirectory = SimulationDirectory;

% Choose hypothetical buoy locations
%arrayShape = 'circle';
arrayShape = 'circle';

switch arrayShape
    case 'circle'
        % Setup buoy array - Circle
        x_center_ind = round(Nx/2);
        y_center_ind = round(Ny/2);
        x_center = x(x_center_ind);
        y_center = y(y_center_ind);
        n_theta = 3;
        theta_spread = 120; % deg
        dir_principal = cell2mat(dspec2char(S_Measured,'Wdir'))*180/pi; % peak direction
        theta = linspace(dir_principal-theta_spread/2,dir_principal+theta_spread/2,n_theta)';
        circle_radius = [100, 150]; % m
        [x_array,y_array,z_array] = setBuoyArray_Circle(x_center,y_center,theta,circle_radius,dx,dy,Y);
        t_array = 1:Nt;
    case 'box'
        % Setup buoy array - box
        dir_principal = cell2mat(dspec2char(S_Measured,'Wdir'))*180/pi; % peak direction
        box_side_length = 200;
        x_center_ind = round(Nx/2);
        y_center_ind = round(Ny/2);
        x_center = x(x_center_ind);
        y_center = y(y_center_ind);
        num_buoys = 4;
        [x_array,y_array,z_array] = setBuoyArray_Box(x_center,y_center,box_side_length,num_buoys,dir_principal,dx,dy,Y);
        t_array = 1:Nt;
end
H_sig = spec2char(S_Measured,'Hm0');


%% Least squares calculation

% set time bursts duration/delay etc
T_meas = 15; %sec
T_pred = 15; %sec
T_delay = 2; %sec
overlap = 0.5;

% set least squares model values
k = logspace(-3,1,51)*2*pi;
theta_wavenumber = linspace(-180,180,19)*pi/180;
theta_wavenumber = theta_wavenumber(2:end);
reg_factor = 5e2;

% Set target point(s)
x_target = (round(Nx/4):1:round(Nx/2))*dx;
y_target = (round(Ny/2-20):round(Ny/2+20))*dy;
%ind = 1;
%x_target = x_array(ind);
%y_target = y_array(ind);

% Calculate least squares prediction
[z_target_pred,z_target_truth,t_pred] = runLeastSquaresPrediction(...
    x_target,y_target,dt,x_array,y_array,z_array,...
    k,theta_wavenumber,reg_factor,T_meas,T_pred,T_delay,overlap,Y);

% [z_target_pred,z_target_truth,t_pred] = runLeastSquaresPrediction_Swifts(...
%     repmat(x_array,[Nt,1]),repmat(y_array,[Nt,1]),z_array',ind,dt,...
%     k,theta_wavenumber,reg_factor,T_meas,T_pred,T_delay,overlap);
%% Compare time series prediction at point
if plotPointComparison
    x_target_ind_plot = 1;%find(x_target == round(Nx/4)*dx);
    y_target_ind_plot = 1;%find(y_target == round(Ny/2)*dy);
    target_ind_plot = (x_target_ind_plot-1)*length(y_target)+y_target_ind_plot;
    t_pred_plot = squeeze(t_pred(:,:,target_ind_plot));
    z_target_pred_plot = squeeze(z_target_pred(:,:,target_ind_plot));
    z_target_truth_plot = squeeze(z_target_truth(:,:,target_ind_plot));
    fig(1) = figure(1);
    clf(fig(1))
    fig(1).Position = [1 1 8 6];
    fig(1).PaperPosition = fig(1).Position;
    subplot(2,1,1)
    hold on
    p1 = plot(t_pred_plot',z_target_truth_plot','-b');
    p2 = plot(t_pred_plot',z_target_pred_plot','-k');
    hold off
    box on
    legend([p1(1),p2(1)],{'Measured','Predicted'})
    ylim([-H_sig H_sig])
    xlim([min([t';t_pred(:)]),max([t';t_pred(:)])])
    xlabel('Time (s)')
    ylabel('\eta (m)')
    
    subplot(2,6,7:8)
    hold on
    pcolor(repmat(Y.x',[Ny,1]),repmat(Y.y,[1,Nx]),Y.Z(:,:,1))
    shading flat
    plot(x_array,y_array,'or')
    plot(x_target(x_target_ind_plot),y_target(y_target_ind_plot),'ok')
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
    print(fig(1),'-dpng',[plotDirectory '/PointPrediction.png'])
    
end
%% Compare surface reconstruction
if plotSurfaceComparison
    surfacePlotDirectory = [plotDirectory '/SurfacePrediction'];
    mkdir(surfacePlotDirectory)
    
    % remove overlap
    [t_unique,ind_unique,~] = unique(t_pred(:,:,1));
    [burst_ind,time_ind] = ind2sub([size(t_pred,1),size(t_pred,2)],ind_unique);
    
    movieInd = find(t_unique>0 & t_unique<60);
    ny_target = length(y_target);
    nx_target = length(x_target);
    for i = 1:length(movieInd)
        fig(2) = figure(2); clf(fig(2));
        fig(2).Position = [2 2 5 8];
        fig(2).PaperPosition = fig(2).Position;
        subplot(2,1,1)
        hold on
        pcolor(x_target,y_target',...
            reshape(squeeze(z_target_truth(burst_ind(movieInd(i)),time_ind(movieInd(i)),:)),[ny_target,nx_target]))       
        shading flat
        plot(x_array,y_array,'or')
        hold off
        set(gca,'CLim',1/2*[-H_sig H_sig],'XLim',[0 Nx*dx],'YLim',[0 Ny*dy])
        xlabel('x [m]'), ylabel('y [m]'),
        cbar = colorbar;
        ylabel(cbar,'\eta (m)')
        box('on')
        title(sprintf('t = %d sec',t_unique(movieInd(i))))
        
        subplot(2,1,2)
        hold on
        pcolor(x_target,y_target,...
            reshape(squeeze(z_target_pred(burst_ind(movieInd(i)),time_ind(movieInd(i)),:)),[ny_target,nx_target]))
        
        shading flat
        plot(x_array,y_array,'or')
        hold off
        set(gca,'CLim',1/2*[-H_sig H_sig],'XLim',[0 Nx*dx],'YLim',[0 Ny*dy])
        xlabel('x [m]'), ylabel('y [m]'),
        cbar = colorbar;
        ylabel(cbar,'\eta (m)')
        box('on')
        print(fig(2),[surfacePlotDirectory '/SurfacePrediction_' sprintf('%03d',i-minInd) '.png'],'-dpng')
    end
end
%% Rendered Frames
if plotRenderedFrames
    RenderedMovieFramesDir = [plotDirectory '/SurfaceMovieFrames_Rendered'];
    mkdir(RenderedMovieFramesDir)
    [xs,ys,zs] = sphere;
    x_target_ind_plot = find(x_target == round(Nx/4)*dx);
    y_target_ind_plot = find(y_target == round(Ny/2)*dy);
    target_ind_plot = (x_target_ind_plot-1)*length(y_target)+y_target_ind_plot;
    t_pred_plot = squeeze(t_pred(:,:,target_ind_plot));
    z_target_pred_plot = squeeze(z_target_pred(:,:,target_ind_plot));
    
    x_target_ind_total = find(Y.x == round(Nx/4)*dx);
    y_target_ind_total = find(Y.y == round(Ny/2)*dy); 
    
    i1 = min(t_pred_plot(:))/dt;
    movieLength = 120; %seconds
    for i = i1:i1+(movieLength/dt)
        fig(3) = figure(3); clf(fig(3));
        fig(3).Position = [2 2 16 10];
        fig(3).PaperPosition = fig(3).Position;
        subplot(2,1,1)
        hold on
        surf(repmat(Y.x',[Ny,1]),repmat(Y.y,[1,Nx]),Y.Z(:,:,i),Y.Z(:,:,i))
        for j = 1:length(x_array)
            ssphere = surf(x_array(j)+xs,y_array(j)+ys,z_array(j,i)+zs,repmat(shiftdim([1 1 0],-1),[21 21 1]));
        end
        surf(Y.x(x_target_ind_total)+xs,Y.y(y_target_ind_total)+ys,squeeze(Y.Z(y_target_ind_total,x_target_ind_total,i))+zs,repmat(shiftdim([1 0 0],-1),[21 21 1]));
        colormap gray
        shading interp
        lighting phong
        material shiny
        hold off
        axis equal
        set(gca,'ZLim',[-H_sig H_sig],'CLim',[-H_sig H_sig],'XLim',...
            [x_target(x_target_ind_plot)-50,x_target(x_target_ind_plot)+250],...
            'YLim',[y_target(y_target_ind_plot)-150,y_target(y_target_ind_plot)+200])
        view(0,15)
        axis off
        light
        subplot(2,1,2)
        hold on
        p1 = plot(Y.t(i1:i),squeeze(Y.Z(y_target_ind_total,x_target_ind_total,i1:i)),'-k');
        p2 = plot(t_pred_plot(t_pred_plot/dt<i),z_target_pred_plot(t_pred_plot/dt<i),'sb');
        box on
        hold off
        legend('Measured','Predicted')
        xlabel('Time (s)')
        ylabel('\eta [m]')
        ylim([-H_sig H_sig]*2/3)
        xlim([i1 i1+120])
        print(fig(3),[RenderedMovieFramesDir '/SimulationFrame_' sprintf('%03d',i-i1+1) '.png'],'-dpng')
    end
end

