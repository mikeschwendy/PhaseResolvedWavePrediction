clc
clearvars
close all

plotDirectory = '/Users/mike/Documents/UW/Research/Results/PapaSwiftVisualization';
swiftFile = '/Users/mike/Documents/UW/Research/Data/StereoTests/SWIFT/allSWIFT2015_withsurfaceelevations.mat';
load(swiftFile)
N = length(allSWIFT);
indPairs = [];
lat = [];
lon = [];
sigwaveheight = [];
peakwaveperiod = [];
peakwavedir = [];
distance = [];
x = [];
y = [];
numPairs = 0;
for ii = 1:N
    %ii
    for jj = (ii+1):N
        if allSWIFT(ii).time == allSWIFT(jj).time
            numPairs = numPairs+1;
            indPairs(numPairs,1) = ii;
            indPairs(numPairs,2) = jj;
            lat(numPairs,1) = allSWIFT(ii).lat;
            lat(numPairs,2) = allSWIFT(jj).lat;
            lon(numPairs,1) = allSWIFT(ii).lon;
            lon(numPairs,2) = allSWIFT(jj).lon;
            sigwaveheight(numPairs,1) = allSWIFT(ii).sigwaveheight;
            sigwaveheight(numPairs,2) = allSWIFT(jj).sigwaveheight;
            peakwaveperiod(numPairs,1) = allSWIFT(ii).peakwaveperiod;
            peakwaveperiod(numPairs,2) = allSWIFT(jj).peakwaveperiod;
            peakwavedir(numPairs,1) = allSWIFT(ii).peakwavedirT;
            peakwavedir(numPairs,2) = allSWIFT(jj).peakwavedirT;
            zone = utmzone([lat(numPairs,1),lon(numPairs,1)]);
            [ellipsoid,estr] = utmgeoid(zone);
            utmstruct = defaultm('utm');
            utmstruct.zone = zone;
            utmstruct.geoid = ellipsoid(1,:);
            utmstruct = defaultm(utmstruct);
            [xii,yii] = mfwdtran(utmstruct,lat(numPairs,1),lon(numPairs,1));
            [xjj,yjj] = mfwdtran(utmstruct,lat(numPairs,2),lon(numPairs,2));
            x(numPairs,1) = xii;
            x(numPairs,2) = xjj;
            y(numPairs,1) = yii;
            y(numPairs,2) = yjj;
            distance(numPairs) = sqrt((xii-xjj).^2+(yii-yjj).^2);
            xc = (xii+xjj)/2;
            yc = (yii+yjj)/2;
            
            if distance(numPairs) < 200 && peakwavedir(numPairs,1) <= 360 &&...
                    peakwavedir(numPairs,1) >= 0 && ~isempty(allSWIFT(ii).z) &&...
                    ~isempty(allSWIFT(jj).z)
                %% Save SWIFT pairs
                SWIFT1 = allSWIFT(ii);
                SWIFT2 = allSWIFT(jj);
                save([plotDirectory '/SWIFTdata_' sprintf('%05d',numPairs) '.mat'],'SWIFT1','SWIFT2')
                
                %% Plot SWIFT data
                
                [~,fPeakIndii] = min(abs(allSWIFT(ii).wavespectra.freq-1/allSWIFT(ii).peakwaveperiod));
                [~,fPeakIndjj] = min(abs(allSWIFT(jj).wavespectra.freq-1/allSWIFT(jj).peakwaveperiod)); 
                fsii = round(length(allSWIFT(ii).z)/(60*9.1));
                fsjj = round(length(allSWIFT(jj).z)/(60*9.1));
                peakWavelength = 9.8/(2*pi)*peakwaveperiod(numPairs,1).^2;
                peakPhasespeed = peakWavelength/peakwaveperiod(numPairs,1);
                averageDelay = distance(numPairs)/peakPhasespeed;
                
                fig(1) = figure(1);
                clf(fig(1))
                subplot(4,1,1)
                plot(x(numPairs,1),y(numPairs,1),'x','markersize',16)
                hold on
                plot(x(numPairs,2),y(numPairs,2),'x','markersize',16)
                plot([xc,xc+distance(numPairs)/2*cos(pi/180*peakwavedir(numPairs))],...
                    [yc,yc+distance(numPairs)/2*sin(pi/180*peakwavedir(numPairs))],'k--')
                plot([xc,xc-distance(numPairs)/2*cos(pi/180*peakwavedir(numPairs))],...
                    [yc,yc-distance(numPairs)/2*sin(pi/180*peakwavedir(numPairs))],'k--')
                title({datestr(allSWIFT(jj).time); ...
                    sprintf('Peak Wavelength = %.1f m',peakWavelength); ...
                    sprintf('Distance = %.1f m',distance(numPairs));...
                    sprintf('Average Delay = %.1f s',averageDelay)})
                hold off
                axis equal
                xlabel('x (m)')
                ylabel('y (m)')
                subplot(4,1,2)
                plot(allSWIFT(ii).wavespectra.freq,allSWIFT(ii).wavespectra.energy)
                hold on
                plot(allSWIFT(jj).wavespectra.freq,allSWIFT(jj).wavespectra.energy)
                plot(allSWIFT(ii).wavespectra.freq(fPeakIndii),allSWIFT(ii).wavespectra.energy(fPeakIndii),'ok')
                plot(allSWIFT(jj).wavespectra.freq(fPeakIndjj),allSWIFT(jj).wavespectra.energy(fPeakIndjj),'ok')
                hold off
                xlabel('Frequency [Hz]')
                ylabel('Energy [m^2/Hz]')
                xlim([0 0.5])
                subplot(4,1,3)
                plot(allSWIFT(ii).wavespectra.freq,180/pi*atan2(allSWIFT(ii).wavespectra.b1,allSWIFT(ii).wavespectra.a1))
                hold on
                plot(allSWIFT(jj).wavespectra.freq,180/pi*atan2(allSWIFT(jj).wavespectra.b1,allSWIFT(jj).wavespectra.a1))
                plot(allSWIFT(ii).wavespectra.freq(fPeakIndii),180/pi*atan2(allSWIFT(ii).wavespectra.b1(fPeakIndii),allSWIFT(ii).wavespectra.a1(fPeakIndii)),'ok')
                plot(allSWIFT(jj).wavespectra.freq(fPeakIndjj),180/pi*atan2(allSWIFT(jj).wavespectra.b1(fPeakIndjj),allSWIFT(jj).wavespectra.a1(fPeakIndjj)),'ok')
                hold off
                xlabel('Frequency [Hz]')
                ylabel('Mean Direction [deg]')  
                ylim([-180 180])
                xlim([0 0.5])
                subplot(4,1,4)
                plot(1/fsii*(1:length(allSWIFT(ii).z)),allSWIFT(ii).z)
                hold on
                plot(1/fsjj*(1:length(allSWIFT(jj).z)),allSWIFT(jj).z)
                hold off
                xlabel('time (s)')
                ylabel('\eta (m)')
                xlim([120 240])
                print(fig(1),[plotDirectory '/SWIFT_DistanceAndWaveDirection_' sprintf('%05d',numPairs) '.jpg'],'-djpeg')
                
                %% Calculate least squares predictions
                t1 = 1/fsii*(1:length(allSWIFT(ii).z));
                z1 = allSWIFT(ii).z;
                nonNan1 = ~isnan(z1);
                z1 = z1(nonNan1);
                t1 = t1(nonNan1);
                t2 = 1/fsii*(1:length(allSWIFT(jj).z));
                z2 = allSWIFT(jj).z;
                nonNan2 = ~isnan(z2);
                z2 = z2(nonNan2);
                t2 = t2(nonNan2);               
                deltaX = distance(numPairs);
                f_oneside = linspace(0.1,0.3,101);
                df = f_oneside(2)-f_oneside(1);
                f = [-f_oneside-df/2, f_oneside];
                %f = [-allSWIFT(jj).wavespectra.freq+.01;allSWIFT(jj).wavespectra.freq];
                [z_pred_2,z_check_1] = leastSquaresWavePropagation(z1,t1,0,t2,deltaX,f);
                [z_pred_1,z_check_2] = leastSquaresWavePropagation(z2,t2,0,t1,deltaX,f);
                ind1 = t1>100 & t1<400;
                ind2 = t2>100 & t2<400;
                
                %% Plot predictions
                fig(2) = figure(2);
                clf(fig(2))
                fig(2).Position = [1 1 12 6];
                fig(2).PaperPosition = fig(2).Position;
                subplot(2,3,1:2)
                plot(t1,z1)
                hold on
                plot(t1,z_pred_1,'-')
                hold off
                xlabel('Time (t)')
                set(gca,'xlim',[100 400],...
                    'ylim',[-sigwaveheight(numPairs,1) sigwaveheight(numPairs,1)])
                legend('Original','Prediction','location','best')

                subplot(2,3,3)
                plot(z1(ind1),z_pred_1(ind1),'.')
                hold on
                plot([-4 4],[-4 4],'--k')
                hold off
                xlabel('Original')
                ylabel('Prediction')
                set(gca,'xlim',[-sigwaveheight(numPairs,1) sigwaveheight(numPairs,1)],...
                    'ylim',[-sigwaveheight(numPairs,1) sigwaveheight(numPairs,1)])
                
                subplot(2,3,4:5)
                plot(t2,z2)
                hold on
                plot(t2,z_pred_2,'-')
                hold off
                xlabel('Time (t)')
                set(gca,'xlim',[100 400],...
                    'ylim',[-sigwaveheight(numPairs,1) sigwaveheight(numPairs,1)])
                legend('Original','Prediction','location','best')

                subplot(2,3,6)
                plot(z2(ind2),z_pred_2(ind2),'.')
                hold on
                plot([-4 4],[-4 4],'--k')
                hold off
                xlabel('Original')
                ylabel('Prediction')
                set(gca,'xlim',[-sigwaveheight(numPairs,1) sigwaveheight(numPairs,1)],...
                    'ylim',[-sigwaveheight(numPairs,1) sigwaveheight(numPairs,1)])
                print(fig(2),[plotDirectory '/SWIFT_LeastSquaresPrediction_' sprintf('%05d',numPairs) '.jpg'],'-djpeg')
                
                %% Test cross-correlations for simple time lag
                [crosscor12,lags12] = xcorr(xcorr(z1,z2));
                [crosscor1,lags1] = xcorr(xcorr(z1(ind1),z_pred_1(ind1)));
                [crosscor2,lags2] = xcorr(xcorr(z2(ind2),z_pred_2(ind2)));
                fig(3) = figure(3);
                clf(fig(3))
                fig(3).Position = [1 1 8 6];
                fig(3).PaperPosition = fig(3).Position;
                subplot(3,1,1)
                plot(lags12*1/fsii,crosscor12)
                xlim([-60 60])
                title('Original 1 / Original 2')
                subplot(3,1,2)
                plot(lags1*1/fsii,crosscor1)
                xlim([-60 60])
                title('Original 1 / Prediction 1')
                subplot(3,1,3)
                plot(lags2/fsjj,crosscor2)                
                xlim([-60 60])
                title('Original 2 / Prediction 2')
                print(fig(3),[plotDirectory '/SWIFT_CrossCorrelation_' sprintf('%05d',numPairs) '.jpg'],'-djpeg')
            end
        end
    end
end