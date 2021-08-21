%% plot smooth TFR data
% ITPC for more than one group!
% data is a matrix of frequency (rows) x time (columns)


function it_plotSmoothTFR(data, timeoi,freqoi,chns,LINES)
%%
TITLES = {'TD: ITPC values for time and frequency';'ASD: ITPC values for time and frequency';...
    'TD-ASD: ITPC values for time and frequency'};

LINES = [0,0.65,1.3,1.95,2.6]; %PARAMS.lines;
N = 500;
%timeReduction = [1:length(timeoi)]; %change according to the NaNs in the data matrix. Try to get rid of as many as possible. 
timeReduction = [100:length(timeoi)-100];
Timeoi = timeoi(timeReduction);
clear Data
%Data = data(:,timeReduction); 

Data = data_bl_ASD(:,timeReduction); %data_bl is baselined data (see VAMP_TFR);
  
    [n, m] = size(Data);
    [x,y] = meshgrid(Timeoi,freqoi); % low-res grid
    [x2,y2] = meshgrid(Timeoi(1):1/N/5:Timeoi(end),freqoi(1):.01:freqoi(end));  %high-res grid
    dataInterp = interp2(x,y,Data, x2,y2, 'linear'); %interpolate up


    figure;
    subplot(4,4,[1 2 2 8])

    f = surf(x2,y2,dataInterp);
    f.EdgeColor = 'none';
    f.FaceColor = 'interp';
    f.FaceLighting = 'gouraud';
    set(gca,'ydir','normal')
    ylabel('Frequency (Hz)')
    xlabel('Time (Sec.)');
    colorbar;
    colormap jet;
    ax = gca;
    %caxis([0.2 0.4]); % VAMP: change color scale -0.5:0.5 for TD, -0.3:0.3 for ASD
%     if group ~= 3
%         caxis([-18 20]);
%     end
    caxis(ax.CLim)
    %caxis([0 40]);
    view(0,90)
    axis tight
    hold on;
    % stimuli lines

    z = get(f,'ZData');
    %set(f,'ZData',z-10);
    z_max = max(max(get(f,'Zdata')));
    hold on;

    for k = 1:length(LINES)
         line([LINES(k), LINES(k)],[y2(1,1),y2(end,1)],[z_max,z_max]...
          ,'Color','w','LineWidth',2,'LineStyle','--');
    end

    set(gca,'fontsize', 14);
    hold on;
    dataInterp(dataInterp==1) = NaN;
    timeVec = Timeoi(1):1/N/5:Timeoi(end);
    ax1 = subplot(4,4,[9, 12]);
    area(ax1,timeVec, nanmean(dataInterp,1),'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
    ylim([0.12 0.5])
    box off
    %plot(timeVec, mean(dataInterp,1),'LineWidth',2);
    %axis tight
    xlabel('Time, sec.')
    ylabel('ITPC value')
    set(gca,'fontsize', 14);
    axis tight


%if plotP ==1
% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(4,4,[1 4 2 8])
%      h =  imagesc(Timeoi(1):1/N/5:Timeoi(end),freqoi(1):.01:freqoi(end),dataInterp)
%     hold on; 
%     dataInterp(isnan(dataInterp)==1)=1;
%     pSig = dataInterp < 0.05;
%     %[row,col,v] = find(pSig);   
%     contour(Timeoi(1):1/N/5:Timeoi(end),freqoi(1):.01:freqoi(end),pSig,'k');
%     colormap jet;
%     colorbar;
%     caxis([-0.1 0.1]);
%     set(gca,'YDir','normal')
%     for k = 2:length(lines)
%         line([lines(k), lines(k)],[y2(1,1),y2(end,1)],[z_max,z_max]...
%         ,'Color','w','LineWidth',2,'LineStyle','--');
%     end
%     ylabel('Frequency (Hz)');
%     xlabel('Time (Sec.)');
%     title('TD-ASD ITPC difference')
%     set(gca,'fontsize', 14);
% 
% hold on; 
% subplot(3,3,[1 4])
% 
% %plot(mean(dataInterp,2),freqoi)
% freqVec = freqoi(1):.01:freqoi(end); 
% %freqVec = freqoi(1):1/500/5:freqoi(end);
% plot(mean(dataInterp,2),freqVec,'LineWidth',2);
% 
% %figure; plot(nanmean(data,2),freqoi,'LineWidth',2);
% 
% %set(gca, 'xdir','reverse')
% axis tight
% ylabel('Frequency, Hz')
% xlabel('ITPC value')
% set(gca,'fontsize', 12);


end