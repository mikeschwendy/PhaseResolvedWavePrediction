clc
clearvars
close all

performPointComparison = true;
plotRenderedFrames = false;
performSurfaceReconstruction = false;

% Load simulation data
SimulationDirectory = '/Users/mike/Dropbox/PhaseResolvedWavePrediction/SurfaceSimulations/CDIP_SanNicolasIsland_Mar2016/CDIP_04Mar2016_L500_dx4_T300_dt1';
load([SimulationDirectory '/SimulationData.mat'],'Y','dx','dt','L','T','Nx','Nt','x','t','Ny','dy','y','S_Measured')
plotDirectory = SimulationDirectory;

% Choose hypothetical buoy locations
arrayShape = 'circle';

switch arrayShape
    case 'circle'
        % Setup buoy array - Circle
        x_center_ind = round(Nx/2);
        y_center_ind = round(Ny/2);
        x_center = x(x_center_ind);
        y_center = y(y_center_ind);
        n_theta = 3;
        theta_spread = 90; % deg
        dir_principal = cell2mat(dspec2char(S_Measured,'Wdir'))*180/pi; % peak direction
        theta = linspace(dir_principal-theta_spread/2,dir_principal+theta_spread/2,n_theta)';
        circle_radius = [150]; % m
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
T_meas = 30; %sec, length of observation window
T_pred = 10; %sec, length of prediction window
T_delay = 5; %sec, 0=prediction window follows immediately after measurement, negative means some prediction happens during measurement (hindcast)
overlap = 0.3; % 1=no overlap, 0.5=50%, etc

% set least squares model values
k = logspace(-3,-1,41)*2*pi;
theta_wavenumber = linspace(-180,180,13)*pi/180;
theta_wavenumber = theta_wavenumber(2:end);
reg_factor = 1e0;

%% Compare time series prediction at point
if performPointComparison
    
    x_target = repmat(x_center,[Nt,1]);
    y_target = repmat(y_center,[Nt,1]);
    z_target = squeeze(Y.Z(y_center_ind,x_center_ind,:));
    
    [z_pred,z_truth,t_pred] = runLeastSquaresPrediction_PointComparison(...
        x_target,y_target,z_target,dt,repmat(x_array',[Nt,1]),repmat(y_array',[Nt,1]),z_array',...
        k,theta_wavenumber,reg_factor,T_meas,T_pred,T_delay,overlap);
    
    
    fig(1) = figure(1);
    clf(fig(1))
    subplot(2,1,1)
    hold on
    p1 = plot(t_pred',z_truth','-b');
    p2 = plot(t_pred',z_pred','-k');
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
    plot(x_target,y_target,'ok')
    hold off
    set(gca,'CLim',[-H_sig H_sig],'XLim',[0 Nx*dx],'YLim',[0 Ny*dy])
    xlabel('x [m]'), ylabel('y [m]'),
    
    subplot(2,4,7:8)
    [z_regress,z_regress_int,z_regress_resid,~,regress_stats] = regress(z_truth(:),z_pred(:));
    regress_r2 = regress_stats(1);
    plot(z_pred(:),z_truth(:),'.')
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
    %% Rendered prediction
    if plotRenderedFrames
        RenderedMovieFramesDir = [plotDirectory '/SurfaceMovieFrames_Rendered'];
        mkdir(RenderedMovieFramesDir)
        [xs,ys,zs] = sphere;
        % remove overlap
        [t_unique,ind_unique,~] = unique(t_pred);
        [burst_ind,time_ind] = ind2sub([size(t_pred)],ind_unique);
        movieInd = find(t_unique>0 & t_unique<180);
        z_truth_plot = z_truth(ind_unique(movieInd));
        z_pred_plot = z_pred(ind_unique(movieInd));
        t_pred_plot = t_pred(ind_unique(movieInd));
        
        for i = 1:length(movieInd)
            fig(3) = figure(3); clf(fig(3));
            fig(3).Units = 'inches';
            fig(3).Position = [2 2 16 10];
            fig(3).PaperPosition = fig(3).Position;
            subplot(2,1,1)
            hold on
            surf(repmat(Y.x',[Ny,1]),repmat(Y.y,[1,Nx]),Y.Z(:,:,movieInd(i)),Y.Z(:,:,movieInd(i)))
            for j = 1:length(x_array)
                ssphere = surf(x_array(j)+xs,y_array(j)+ys,z_array(j,movieInd(i))+zs,repmat(shiftdim([1 1 0],-1),[21 21 1]));
            end
            surf(x_target(movieInd(i))+xs,y_target(movieInd(i))+ys,z_target(movieInd(i))+zs,repmat(shiftdim([1 0 0],-1),[21 21 1]));
            colormap gray
            shading interp
            lighting phong
            material shiny
            hold off
            axis equal
            set(gca,'ZLim',[-H_sig H_sig],'CLim',[-H_sig H_sig],...
                'XLim',[x_target(1)-50,x_target(1)+250],...
                'YLim',[y_target(1)-150,y_target(1)+200])
            view(0,15)
            axis off
            light
            subplot(2,1,2)
            hold on
            p1 = plot(t_pred_plot(1:i),z_truth_plot(1:i),'-');
            p2 = plot(t_pred_plot(1:i),z_pred_plot(1:i),'-');
            box on
            hold off
            legend('Measured','Predicted')
            xlabel('Time (s)')
            ylabel('\eta [m]')
            ylim([-H_sig H_sig]*2/3)
            xlim([min(t_pred_plot) max(t_pred_plot)])
            print(fig(3),[RenderedMovieFramesDir '/SimulationFrame_' sprintf('%03d',i) '.png'],'-dpng')
        end
    end
end
%% Make surface reconstruction
if performSurfaceReconstruction
    x_target = (round(Nx/4):1:round(Nx/2))*dx;
    y_target = (round(Ny/2-20):round(Ny/2+20))*dy;
    [x_target,y_target] = meshgrid(x_target,y_target);
    
    [z_target_pred,t_pred] = runLeastSquaresPrediction_SurfaceReconstruction(...
        x_target,y_target,repmat(x_array',[Nt,1]),repmat(y_array',[Nt,1]),z_array',dt,...
        k,theta_wavenumber,reg_factor,T_meas,T_pred,T_delay,overlap);
    
    surfacePlotDirectory = [plotDirectory '/SurfacePrediction'];
    mkdir(surfacePlotDirectory)
    
    % remove overlap
    [t_unique,ind_unique,~] = unique(t_pred(:,:,1));
    [burst_ind,time_ind] = ind2sub([size(t_pred,1),size(t_pred,2)],ind_unique);
    
    movieInd = find(t_unique>0 & t_unique<180);
    [ny_target,nx_target] = size(x_target);
    for i = 1:length(movieInd)
        fig(2) = figure(2); clf(fig(2));
        fig(2).Units = 'inches';
        fig(2).Position = [2 2 5 8];
        fig(2).PaperPosition = fig(2).Position;
        subplot(2,1,1)
        hold on
        pcolor(Y.x,Y.y',Y.Z(:,:,t_unique(movieInd(i))))
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
        print(fig(2),[surfacePlotDirectory '/SurfacePrediction_' sprintf('%03d',i) '.png'],'-dpng')
    end
end
