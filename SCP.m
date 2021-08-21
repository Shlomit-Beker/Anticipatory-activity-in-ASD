% datamatrix1 and datamatrix2 in this case are [subjects x time x electrodes] for each condition
% Shlomit Beker 2017 (following a version from the CNL)

function ST = SCP(datamatrix,newChansInds, figurenum, time,a,sigVals)

datamatrix1 = datamatrix{1};
datamatrix2 = datamatrix{2};
% Initial parameters
alpha_level = 0.05;
numconsectime = 2; %number of consecutive time points required
numconsecchan = 2; %number of consecutive channels required
No_of_Electrodes = 64; %number of electrodes


% Run through all the t-tests
for q = 1:length(datamatrix1(1,:,1)) %this loops through each time point
    for r = 1:length(datamatrix1(1,1,:)) % this loops through each electrode
       [H(q,r) P(q,r) C D(q,r)]=ttest2(datamatrix1(:,q,r),datamatrix2(:,q,r));       
       t_val(q,r) = D(q,r).tstat;
    end
end

% Threshold t-stat values
[P_val T_val]=SCP_THRESHOLD(squeeze(P)',t_val',alpha_level,numconsectime,numconsecchan);


% Reorder the channels to an anatomical order
newT_val = T_val(newChansInds,:);
newsigVals = sigVals(newChansInds,:);

% Make the figure
sigImage = zeros(64,308);
sigImage(newsigVals) = newT_val(newsigVals);
figure;
imagesc(time,[1:No_of_Electrodes],sigImage);
% cmap = jet;
% cmap(end,:) = 1;
caxis([-5 5]);
colormap(bluewhitered(256)), colorbar
hold on;
imagesc(time,[1:No_of_Electrodes],squeeze(newT_val),'AlphaData',0.3);
colormap(bluewhitered(256)), colorbar

caxis([-5 5]);
colorbar;
ylabel('Electrode index')
xlabel('Time (ms)');
y1=get(gca,'ylim');
mini = min(abs(time));
hold on;
plot([mini mini],y1,'--k');
set(gca,'fontsize', 12);



