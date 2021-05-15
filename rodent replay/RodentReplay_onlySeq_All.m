%% Explicit model time in the first level GLM
close all
clear;
clc;

%% prepare data
load SpikingSleep
nshuf     = 100;
nlags     = 600; % in 1 ms lag
binSize   = .001; % 1 ms

%% TDLM in multi-scale (time lag of interest, 2 ms -> 10m/s) 
X=[];
for ievent=1:length(Events) % loop over all events
    X=[X;Events{ievent}]; % ntime (concatenated) * ncells
end

[Seq_integrate, Seq_separate]=TDLM_multiscale(X, RateMap, binSize, nshuf, nlags);
save('OnlySeq_full_weighted','Seq_integrate', 'Seq_separate', '-v7.3')

%% Stats & Plot
load OnlySeq_full_weighted

%% sequence effect
figure
subplot(2,2,1)
xx= squeeze(Seq_integrate(1,1,:,:));
plot(xx(1,:))
hold on
yline(prctile(max(squeeze(xx(2:end,:)),[],2),95),'--r');
title('In-bound: Fwd Seq')
xlabel(' time lag (1 ms)')
ylabel('sequence strength')
ylim([-0.002,0.007])

subplot(2,2,2)
xx= squeeze(Seq_integrate(1,2,:,:));
plot(xx(1,:))
hold on
yline(prctile(max(squeeze(xx(2:end,:)),[],2),95),'--r');
title('Out-bound: Fwd Seq')
xlabel(' time lag (1 ms)')
ylabel('sequence strength')
ylim([-0.002,0.007])

subplot(2,2,3)
xx= squeeze(Seq_integrate(2,1,:,:));
plot(xx(1,:))
hold on
yline(prctile(max(squeeze(xx(2:end,:)),[],2),95),'--r');
title('In-bound: Bkw Seq')
xlabel(' time lag (1 ms)')
ylabel('sequence strength')
ylim([-0.002,0.007])

subplot(2,2,4)
xx= squeeze(Seq_integrate(2,2,:,:));
plot(xx(1,:))
hold on
yline(prctile(max(squeeze(xx(2:end,:)),[],2),95),'--r');
title('Out-bound: Bkw Seq')
xlabel(' time lag (1 ms)')
ylabel('sequence strength')
ylim([-0.002,0.007])

%% sequence effect -> fwd+bkw, fwd-bkw
figure
subplot(1,2,1)
xx= squeeze(Seq_integrate(1,1,:,:)) + squeeze(Seq_integrate(2,1,:,:));
plot(xx(1,:))
hold on
yline(prctile(max(squeeze(xx(2:end,:)),[],2),95),'--r');
title('In-bound: Fwd + Bkw Seq')
xlabel(' time lag (1 ms)')
ylabel('sequence strength')
ylim([-0.002,0.01])

subplot(1,2,2)
xx= squeeze(Seq_integrate(1,2,:,:)) + squeeze(Seq_integrate(2,2,:,:));
plot(xx(1,:))
hold on
yline(prctile(max(squeeze(xx(2:end,:)),[],2),95),'--r');
title('Out-bound: Fwd + Bkw Seq')
xlabel(' time lag (1 ms)')
ylabel('sequence strength')
ylim([-0.002,0.01])

%% sequence effect -> fwd-bkw
figure
subplot(1,2,1)
xx= squeeze(Seq_integrate(1,1,:,:)) - squeeze(Seq_integrate(2,1,:,:));
plot(xx(1,:))
hold on
yline(prctile(max(abs(squeeze(xx(2:end,:))),[],2),95),'--r');
yline(prctile(max(abs(squeeze(xx(2:end,:))),[],2),95)*-1,'--r');

title('In-bound: Fwd - Bkw Seq')
xlabel(' time lag (1 ms)')
ylabel('sequence strength')
ylim([-0.006,0.006])

subplot(1,2,2)
xx= squeeze(Seq_integrate(1,2,:,:)) - squeeze(Seq_integrate(2,2,:,:));
plot(xx(1,:))
hold on
yline(prctile(max(abs(squeeze(xx(2:end,:))),[],2),95),'--r');
yline(prctile(max(abs(squeeze(xx(2:end,:))),[],2),95)*-1,'--r');

title('Out-bound: Fwd - Bkw Seq')
xlabel(' time lag (1 ms)')
ylabel('sequence strength')
ylim([-0.006,0.006])
