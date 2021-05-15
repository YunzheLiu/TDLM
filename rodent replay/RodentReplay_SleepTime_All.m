%% Explicit model time in the first level GLM
close all
clear;
clc;

%% prepare data
load SpikingSleep
nshuf     = 100;
binSize   = .001; % 1 ms
TOI       = 2; % time lag of interest, e.g., 2 ms -> 10m/s
nlags     = TOI*5; % in 1 ms lag

%% TDLM in multi-scale (time lag of interest) 
X=[];
for ievent=1:length(Events) % loop over all events
    X=[X;Events{ievent}]; % ntime (concatenated) * ncells
end
[Seq_integrate, Seq_separate, Time_integrate, Time_separate]=TDLM_multiscale_modelT(X, RateMap, binSize, nshuf, nlags,TOI);
save(['TimeSeq_',num2str(TOI),'ms'],'Seq_integrate', 'Seq_separate', 'Time_integrate', 'Time_separate')

%% Stats & Plot
load (['TimeSeq_',num2str(TOI),'ms'])

%% sequence effect
figure
subplot(2,2,1)
xx= squeeze(Seq_integrate(1,1,:));
histogram(xx(2:end),10,'Normalization','probability')
hold on
xline(prctile(xx(2:end),97.5),'--k');
xline(prctile(xx(2:end),2.5),'--k');
xline(xx(1),'-r');
title(['In-bound: Fwd Seq (', num2str(round(0.02/TOI*1000,2)), ' m/s)'])

subplot(2,2,2)
xx= squeeze(Seq_integrate(1,2,:));
histogram(xx(2:end),10,'Normalization','probability')
hold on
xline(prctile(xx(2:end),97.5),'--k');
xline(prctile(xx(2:end),2.5),'--k');
xline(xx(1),'-r');
title(['Out-bound: Fwd Seq (', num2str(round(0.02/TOI*1000,2)), ' m/s)'])

subplot(2,2,3)
xx= squeeze(Seq_integrate(2,1,:));
histogram(xx(2:end),10,'Normalization','probability')
hold on
xline(prctile(xx(2:end),97.5),'--k');
xline(prctile(xx(2:end),2.5),'--k');
xline(xx(1),'-r');
title(['In-bound: Bkw Seq (', num2str(round(0.02/TOI*1000,2)), ' m/s)'])

subplot(2,2,4)
xx= squeeze(Seq_integrate(2,2,:));
histogram(xx(2:end),10,'Normalization','probability')
hold on
xline(prctile(xx(2:end),97.5),'--k');
xline(prctile(xx(2:end),2.5),'--k');
xline(xx(1),'-r');
title(['Out-bound: Bkw Seq (', num2str(round(0.02/TOI*1000,2)), ' m/s)'])

suptitle('sequence effect')

%% time effect
figure
subplot(2,2,1)
xx= squeeze(Time_integrate(1,1,:));
histogram(xx(2:end),10,'Normalization','probability')
hold on
xline(prctile(xx(2:end),97.5),'--k');
xline(prctile(xx(2:end),2.5),'--k');
xline(xx(1),'-r');
title(['In-bound: Fwd Seq (', num2str(round(0.02/TOI*1000,2)), ' m/s)'])

subplot(2,2,2)
xx= squeeze(Time_integrate(1,2,:));
histogram(xx(2:end),10,'Normalization','probability')
hold on
xline(prctile(xx(2:end),97.5),'--k');
xline(prctile(xx(2:end),2.5),'--k');
xline(xx(1),'-r');
title(['Out-bound: Fwd Seq (', num2str(round(0.02/TOI*1000,2)), ' m/s)'])

subplot(2,2,3)
xx= squeeze(Time_integrate(2,1,:));
histogram(xx(2:end),10,'Normalization','probability')
hold on
xline(prctile(xx(2:end),97.5),'--k');
xline(prctile(xx(2:end),2.5),'--k');
xline(xx(1),'-r');
title(['In-bound: Bkw Seq (', num2str(round(0.02/TOI*1000,2)), ' m/s)'])

subplot(2,2,4)
xx= squeeze(Time_integrate(2,2,:));
histogram(xx(2:end),10,'Normalization','probability')
hold on
xline(prctile(xx(2:end),97.5),'--k');
xline(prctile(xx(2:end),2.5),'--k');
xline(xx(1),'-r');
title(['Out-bound: Bkw Seq (', num2str(round(0.02/TOI*1000,2)), ' m/s)'])

suptitle('time effect')

%% comparing fwd-bkw effect, 
figure

subplot(1,2,1)
xx= squeeze(Seq_integrate(1,2,:)) - squeeze(Seq_integrate(2,2,:));
histogram(xx(2:end),10,'Normalization','probability')
hold on
xline(prctile(xx(2:end),97.5),'--k');
xline(prctile(xx(2:end),2.5),'--k');
xline(xx(1),'-r');
title(['Sequence effect (Out-bound): fwd vs. bkw (', num2str(round(0.02/TOI*1000,2)), ' m/s)'])

subplot(1,2,2)
xx= squeeze(Time_integrate(1,2,:)) - squeeze(Time_integrate(2,2,:));
histogram(xx(2:end),10,'Normalization','probability')
hold on
xline(prctile(xx(2:end),97.5),'--k');
xline(prctile(xx(2:end),2.5),'--k');
xline(xx(1),'-r');
title(['Time effect (Out-bound): fwd vs. bkw (', num2str(round(0.02/TOI*1000,2)), ' m/s)'])


