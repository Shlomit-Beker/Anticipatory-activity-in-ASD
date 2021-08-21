
% Phase Locking (Intra-Trial-Phase-Coherence) analysis. Shlomit Beker 2018
% runs on selected channels. Run this for correlation with behavior for
% Anticipation-in-autism analysis. 
% Shlomit Beker 2019-2020

%% parameters for phase locking
clear PLallstim PLV
CHAN = {'C1','C2','Cz'}; % for entrainment: 'C1','Cz','C2' for visual sequence: 'AF3','Fp1','AFz','Fpz','A7','PO3','POz','PO4',
C = [find(strcmp(ERPb{1}{1}.label,CHAN{1})),find(strcmp(ERPb{1}{1}.label,CHAN{2})),find(strcmp(ERPb{1}{1}.label,CHAN{3}))];%...
CHANNELS = C; % **PLV: C channels; Phase: POz,PO3,PO4

SAMP_RATE = 256; % of the trials, after preprocessing
LOW_FREQUENCY = 0.6;
HIGH_FREQUENCY = 2.5; %range of frequencies on which to make the coherence
TITLES = {'TD Cue','TD No cue','ASD Cue','ASD No Cue'};
OMEGA = 6;
LENGTH_WIN = 3.5; %length of trial in sec
TIME_WINDOW = 1:LENGTH_WIN*SAMP_RATE;
COLORS = {'k',[122,122,122]./255,[212,32,143]./255,[170,125,154]./255};
RAND = 1;
clear i;
FOI = 18; %location of 1.5Hz in the frequencies vector
start = 400;
STIM_TIMES = [start,start+650,start+650*2,start+650*3]/1000*256;

prompt_plotPhase = 'Plot Phases ? (0-no, 1-ind stim, 2-across stim) ';
plotPhase = input(prompt_plotPhase);
if plotPhase == 1
    prompt_stimNum = 'Pick stim (1-4) ';
    stimPick = input(prompt_stimNum);
end

%% Create the data mat files, by condition
tic
flag = 0;
clear sumAngles STphase1 Angles Phase
%PLmat = cell(1,10); % the final cell array will include the four conditions (cue*group)
Angles = [];
participantFlag = 0;
for COND = 1:length(DATA)
    for stim = 1:length(STIM_TIMES)
        for participant = 1:length(DATA{COND})
            clear STphase1
            currentData = DATA{COND}{participant};
            cond = currentData;
            Ntrial = length(cond);
            for i_trial = 1:Ntrial
                if size(cond{i_trial},2) == size(cond{1},2) %to control shorter trials (ignore them)
                    flag = flag+1;
                    [wave,period,scale,cone_of_influence] = basewave4(squeeze(mean(cond{i_trial}(CHANNELS,:)))',SAMP_RATE,LOW_FREQUENCY,HIGH_FREQUENCY,OMEGA,0);
                    STphase1(:,:,i_trial) = squeeze(angle(wave));
                    %STphase1_degrees(:,:,i_trial) = radtodeg(STphase1(:,:,i_trial));%convert to degrees
                end
            end
            frequencies = 1./period;
            
            TOI = ceil(STIM_TIMES(stim));
            %190:192 for 2.5 epoch; %times around peaks (4th stimuli time)
            
            PL = squeeze(mean((abs(sum(exp(1i*STphase1(:,TIME_WINDOW,:)),3))/Ntrial),2));
            PLallstim(:) = PL;
            %%%%%%%%%%%%%%%%%%%%%
%             count = 0;
%             for rand_stim = 1:5
 %                randPL(rand_stim,:) = squeeze(mean((abs(sum(exp(1i*STphase_all_stim(:,TIME_WINDOW,count+[1:Ntrial(rand_stim)])),3))/Ntrial(rand_stim)),2));
%                 count = Ntrial(rand_stim);
%             end
            
            Phase{stim,COND}{participant} = squeeze(STphase1(FOI,TOI,:))';
            PLV{COND}{participant} = PLallstim;
            
        end
        %meanSumAngles(COND) = mean(sumAngles{COND});
        
    end
end
%% plot phases for individual participants, in a selected stimulus location (stimPick)
TITLES = {'TD Cue','TD No cue','ASD Cue','ASD No Cue'};
if plotPhase ~=0 
    for ii = 1:length(DATA)
        figure(ii);
        %suptitle(TITLES{ii});
        hold on;
        for jj = 1:length(DATA{ii}) %participant
            subplot(5,8,jj)
            polarhistogram(Phase{stimPick,ii}{jj},18,'FaceColor',COLORS{ii});
            %polarhistogram(Phase{stim,COND}{participant},18,'FaceColor',COLORS{COND});
            %pause
            title(num2str(ERPb{ii}{jj}.name))
            hold on;
        end
    end
end
%% polar hists of phases for each stimulus timing
   
for ii = 1:size(Phase,1)
    Fig_s = figure(100+ii); title([num2str(ii)]);
    for jj = 1:size(Phase,2)
        for kk = 1:length(Phase{ii,jj})
            Phase_deg{ii,jj}(kk) = meanangle(rad2deg(Phase{ii,jj}{kk})); %average phase, per participant, across trials
            Phase_rad{ii,jj}(kk) = deg2rad(Phase_deg{ii,jj}(kk)); % the same as above, in rads.
        end
         figure(Fig_s); subplot(2,2,jj);
         polarhistogram(Phase_rad{ii,jj},10,'FaceColor',COLORS{jj});
         hold on; 
         %Phase_allStim_trials{jj} = cat(1,Phase{1,jj},Phase{2,jj},Phase{3,jj},Phase{4,jj});
    end
end

%% phases across stim, per participant, plus plot
clear Phase_allStimSubject Mat_1 Mat_2
for ii = 1:size(Phase, 2)
   Fig_s = figure(10+ii); title([num2str(ii)]);

   Mat_1{ii} = cat(1,Phase{1,ii},Phase{2,ii},Phase{3,ii},Phase{4,ii});
     for jj = 1:size(Mat_1{ii},2)
        Mat_2{ii}{jj} = cell2mat(Mat_1{ii}(1:size(Phase, 1),jj));
        Phase_allStimSubject{ii}{jj} = meanangle(rad2deg(Mat_2{ii}{jj}),1);
         figure(Fig_s); subplot(5,7,jj);
         polarhistogram(Phase_allStimSubject{ii}{jj},18,'FaceColor',COLORS{ii});
         hold on; 
     end
end

%%
clear Phase_deg_sum Phase_deg_sum_rad
Fig_a = figure(200); suptitle('all cues');
Fig_m = figure(300); suptitle('mean cue');
for ii=1:4
        Phase_rad_sum{ii} = cat(2,Phase_rad{1,ii},Phase_rad{2,ii},Phase_rad{3,ii},Phase_rad{4,ii});
        figure(Fig_a); subplot(2,2,ii); polarhistogram(Phase_rad_sum{ii},10,'FaceColor',COLORS{ii});title(TITLES{ii});
        hold on; 
        %polarplot([0 real(zm(ii))], [0, imag(zm(ii))],'r');
        hold on;
        Phase_deg_sum{ii} = meanangle(cat(1,Phase_deg{1,ii},Phase_deg{2,ii},Phase_deg{3,ii},Phase_deg{4,ii}));
        %Phase_deg_sum{ii} = meanangle(cat(1,Phase_deg{2,ii},Phase_deg{3,ii},Phase_deg{4,ii}));

        Phase_deg_sum_rad{ii} = deg2rad(Phase_deg_sum{ii});
        figure(Fig_m); subplot(2,2,ii); polarhistogram(Phase_deg_sum_rad{ii},20,'FaceColor',COLORS{ii});
        hold on; 
end

%% Plotting circular histograms using circ statistics toolbox. By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

figure;
subplot(2,2,1)
[a, phi(1),zm(1)] = circ_plot(Phase_deg_sum_rad{1}','hist',[],20,true,true,'linewidth',2,'color','r'); title(TITLES{1}); set(gca,'fontsize', 14);
subplot(2,2,2)
[a, phi(2),zm(2)] = circ_plot(Phase_deg_sum_rad{2}','hist',[],20,true,true,'linewidth',2,'color','r'); title(TITLES{2}); set(gca,'fontsize', 14);
subplot(2,2,3)
[a, phi(3),zm(3)] = circ_plot(Phase_deg_sum_rad{3}','hist',[],20,true,true,'linewidth',3,'color','r'); title(TITLES{3}); set(gca,'fontsize', 14);
subplot(2,2,4)
[a, phi(4),zm(4)] = circ_plot(Phase_deg_sum_rad{4}','hist',[],20,true,true,'linewidth',4,'color','r'); set(gca,'fontsize', 14); title(TITLES{4}); 



%% calculating mean resultant vector length and Rayleigh test
for ii = 1:length(Phase_deg_sum_rad)
    R(ii) = circ_r(Phase_deg_sum_rad{ii}');
    p_alpha(ii) = circ_rtest(Phase_deg_sum_rad{ii}')
end


%% calculate the vector length for each participant, to make anova test
R_group=[];
for ii = 1:length(Phase_allStimSubject)
    for jj = 1:length(Phase_allStimSubject{ii})
       R_group{ii}(jj) = circ_r(deg2rad(Phase_allStimSubject{ii}{jj}'));
    end
end
%% Kuiper test for cdf of samples distributions
[pval, k, K] = circ_kuipertest(Phase_deg_sum_rad{1}, Phase_deg_sum_rad{3}, 100, 1)

%% Watson-Williams test for equality of mean directions (relevant only if the phase itself is important). 
[pval table] = circ_wwtest(Phase_deg_sum_rad{1}, Phase_deg_sum_rad{3})

%% resultant vector length should be > 0.7
[pval, f] = circ_ktest(Phase_deg_sum_rad{1}, Phase_deg_sum_rad{3})

%% [pval, stats] = circ_hktest(alpha, idp, idq, inter, fn)

alpha1 = Phase_deg_sum_rad{1};
alpha2 = Phase_deg_sum_rad{3};
[pval, stats] = circ_hktest(cat(2,alpha1,alpha2),...
    [1:length(alpha1)], [length(alpha1)+1:length(alpha1)+length(alpha2)], 1, {'TD','ASD'})
%%
[rho pval] = circ_corrcc(alpha1, alpha2)

%% make all phases positive (adding 360 for deg / 2*pi for rad). 
for ii = 1:length(Phase_deg_sum_rad)
    for jj = 1:length(Phase_deg_sum_rad{ii})
       if Phase_deg_sum_rad{ii}(jj)<0
           Phase_deg_sum_rad{ii}(jj) = Phase_deg_sum_rad{ii}(jj)+2*pi;
       end
    end
end

% bar graphs with phase per participant. 

figure; 
for ii = 1:length(Phase_deg_sum_rad)
    subplot(2,2,ii);
    bar(Phase_deg_sum_rad{ii},'FaceColor',COLORS{ii});
    title(TITLES{ii});
    hold on; 
    ylim([0 8]);
    ax = gca;
    %ax.YDir = 'reverse'

        xlim([0 1+length(Phase_deg_sum_rad{ii})]);
        hold on;
        xlabel('subject #');
        ylabel('\theta (rad)');
        set(gca,'FontSize',16);
        
        means(ii) = circ_mean(Phase_deg_sum_rad{ii}');
        text(6,7,['\theta = ',num2str(means(ii)+2*pi),'rad',...
            ' (',num2str(round(rad2deg(means(ii)+2*pi))),'\circ',') ']);
end

[p,h,stats] = ranksum(Phase_deg_sum_rad{1},Phase_deg_sum_rad{3})
[p,h,stats] = ranksum(Phase_deg_sum_rad{2},Phase_deg_sum_rad{4})


%% plot (no SD)
figure; 
visTd = cell2mat(PLV{1}');
plot(frequencies, mean(visTd,1),'Color',COLORS{1},'LineWidth',2.5);
novisTD = cell2mat(PLV{2}');
hold on; plot(frequencies, mean(novisTD,1),'Color',COLORS{2},'LineWidth',2.5);
visASD = cell2mat(PLV{3}');
hold on; plot(frequencies, mean(visASD,1),'Color',COLORS{3},'LineWidth',2.5);
novisASD = cell2mat(PLV{4}');
hold on; plot(frequencies, mean(novisASD,1),'Color',COLORS{4},'LineWidth',2.5);
xlabel('Frequency (Hz)');
ylabel('Coherence (AU)');
%xlim([1 2.5]);
%ylim([0.1 0.3]);

legend('TD Cue','TD No Cue','ASD Cue','ASD No Cue');

title('Phase locking values for all conditions')

set(gca,'fontsize', 14);

%% plot individual subjects PLV
%looks better on CP1,CP2,CPz or C1,Cz,C2
figure;
conds = [1,3];
for ii = 1:length(conds)
    for jj = 1:length(PLV{conds(ii)})
        subplot(1,2,ii);
        plot(frequencies, PLV{conds(ii)}{jj},'LineWidth',1.5);
        hold on;
    end
    xlim([1 2.5]);
    ylim([0 3.5]);
    set(gca,'fontsize', 14);
    xlabel('Frequency (Hz)');
    ylabel('Coherence (AU)');
end

%% plot PLV with bounded lines
%looks better on CP1,CP2,CPz or C1,C2,Cz
figure;
for ii = 1:length(PLV)
    %for jj = 1:length(PLV{ii})
    [X,mean_smooth, error_smooth] = drawBoundedLines_NEW(mean(cell2mat(PLV{ii}')),std(cell2mat(PLV{ii}'))./sqrt(length(PLV{ii})),FS,x); % Draw bounded lines
    %[X,mean_smooth, error_smooth] = drawBoundedLines_NEW(mean(cell2mat(PLV{ii}')),std(cell2mat(PLV{ii}'))./sqrt(length(PLV{ii})),FS,x); % Draw bounded lines

    shadedErrorBar(frequencies,mean_smooth,error_smooth,{COLORS{ii},'LineWidth',1},1);
    %pause
    hold on;
   
end

%legend('TD Cue','TD No Cue','ASD Cue','ASD No Cue');
xlabel('Frequency (Hz)');
ylabel('Coherence (AU)');
xlim([LOW_FREQUENCY HIGH_FREQUENCY]);

set(gca,'fontsize', 14);
title('Median Phase locking values for all conditions')

%rank test to test statistical significance in the PLVs
x = mean(cell2mat(PLV{1}'));
y = mean(cell2mat(PLV{3}'));
[p,h,stats] = ranksum(x,y)

% permutation test on two groups
for m = 1:4
    for n = 1:length(PLV{m})
        PLVmeans{m}(n) = mean(PLV{m}{n});
        PLVmax{m}(n) = max(PLV{m}{n}(16:22));
    end
end

%permutation test between groups (cue only)
[p, observeddifference, effectsize] = permutationTest(PLVmeans{1}, PLVmeans{3}, 10000,  'sidedness','larger','plotresult',1)

%max value
[p, observeddifference, effectsize] = permutationTest(PLVmax{1}, PLVmax{3}, 10000,  'sidedness','larger','plotresult',1)

%2-way anova interaction  betweem max value at the area of the ITPC peak

Y = cat(2,PLVmax{1}, PLVmax{2},PLVmax{3}, PLVmax{4}); 
g1 = [ones(1,length(PLVmax{1})),ones(1,length(PLVmax{2})),ones(1,length(PLVmax{3}))*2,ones(1,length(PLVmax{4}))*2];%group 1-TD 2-ASD
g2 =  [ones(1,length(PLVmax{1})),ones(1,length(PLVmax{2}))*2,ones(1,length(PLVmax{3})),ones(1,length(PLVmax{4}))*2];%condition 1-cue 2-noCue
[~,~,stats]= anovan(Y,{g1 g2},'model','interaction','varnames',{'g1','g2'})
results = multcompare(stats,'Dimension',[1 2])
%% plot the PLVmax as a scatter plot (with Mick's function) - requested by Reviewer 3 in round 2 JNP. 
st_boxdotplot([1:4],PLVmax,COLORS([1,4,5,8],:),'iqr',[],[],[],0.3,35,0.5,[],1);

%% calculate differences betweeen a group and zero 
[p, observeddifference, effectsize] = permutationTest(PLVmax{4}, zeros(size(PLVmax{4})), 10000,  'sidedness','larger','plotresult',1)


