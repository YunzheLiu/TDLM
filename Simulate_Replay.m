clear;
clc;
close all;
rng('shuffle')

%% Spects
TF = [0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0]; % transition matrix
TR = TF';
nSensors = 273;
nTrainPerStim = 18; % how many training examples for each stimulus
nNullExamples = nTrainPerStim*8; % how many null examples to use
nSamples = 6000; % 60 seconds of unlabelled data to predict
nSequences = 2000; % how many real sequences to put in the data
maxLag = 60; % evaluate time lag up to 600ms
cTime = 0:10:maxLag*10; % the milliseconds of each cross-correlation time lag
nSubj = 24; % number of subjects to simulate
gamA = 10;
gamB = 0.5;% parameters for the gamma distribution of intervals between states in a sequence
[~, pInds] = uperms([1:8],29);
uniquePerms=pInds;
nShuf = size(uniquePerms,1);
samplerate=100;
nstates=8;

sf = cell(nSubj,1);  sb = cell(nSubj,1);
sf2 = cell(nSubj,1);  sb2 = cell(nSubj,1);

%% Core function
parfor iSj = 1:nSubj
    sf{iSj} = nan(1, nShuf, maxLag+1);
    sb{iSj} = nan(1, nShuf, maxLag+1);
      
    sf2{iSj} = nan(1, nShuf, maxLag+1);
    sb2{iSj} = nan(1, nShuf, maxLag+1);
    
    disp(['iSj=' num2str(iSj)])
    
   %% generate dependence of the sensors
    A = randn(nSensors);
    [U,~] = eig((A+A')/2); 
    covMat = U*diag(abs(randn(nSensors,1)))*U';
    
    %% generate the true patterns
    commonPattern = randn(1,nSensors);    
    patterns = repmat(commonPattern, [8 1]) + randn(8, nSensors);  

    %% make training data
    trainingData = 4*randn(nNullExamples+8*nTrainPerStim, nSensors) + [zeros(nNullExamples,nSensors); ...
        repmat(patterns, [nTrainPerStim 1])];
    trainingLabels = [zeros(nNullExamples,1); repmat((1:8)', [nTrainPerStim 1])];    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Add more noise to some of the patterns
    MoreNoiseind=randsample([1:8],4);
    indend   = MoreNoiseind*nTrainPerStim;
    indstart = (MoreNoiseind-1)*nTrainPerStim+1;
    
    nindex=[];
    
    for iind=1:length(MoreNoiseind)
        xtemp=indstart(iind):indend(iind);
        nindex=[nindex,xtemp];
    end
    
    trainingData(nNullExamples+nindex,:)=trainingData(nNullExamples+nindex,:)+ randn(length(nindex), nSensors);

    %% train classifiers on training data
    betas = nan(nSensors, 8); intercepts = nan(1,8);
    
    for iC=1:8
        [betas(:,iC), fitInfo] = lassoglm(trainingData, trainingLabels==iC, 'binomial', 'Alpha', 1, 'Lambda', 0.006, 'Standardize', false);
        intercepts(iC) = fitInfo.Intercept;
    end
    
    %% make long unlabelled data with or without real sequences in it
    X = nan(nSamples, nSensors);
    X(1,:) = randn([1 nSensors]);
    
    for iT=2:nSamples
         X(iT,:) = 0.95*(X(iT-1,:) + mvnrnd(zeros(1,nSensors), covMat));% add dependence of the sensors 
    end
        
    %% Injecting Sequence
    for iRS = 1:nSequences
        seqTime = randi([40 nSamples-40]); % pick a random point, not right at the edges
%         seqTime = randsample(peakindex,1); % Initiate each sequence at peak of alpha        
        state = false(8,1); state(randi(8)) = true;  % start the sequence in a random state
        
        for iMv=1:2
            if sum(state)==0
                X(seqTime,:) = X(seqTime,:);
            else
                X(seqTime,:) = X(seqTime,:) + patterns(state,:); 
                state = (state'*TR)'; state2 = false(8,1); state2(find(rand < cumsum(state), 1, 'first')) = true; state = state2; % advance states
                seqTime = seqTime + round(gamrnd(gamA,gamB));
            end
        end
    end
    
    %% make predictions with trained models
    preds = 1./(1+exp(-(X*betas + repmat(intercepts, [nSamples 1]))));
    
    %% make decoded time course autocorrelated;      
    coreP = nan(nSamples, 1);    
    coreP(1,1) = randn(1);    

    for iT=2:nSamples
         coreP(iT,1) = 0.95*(coreP(iT-1,1) + randn(1));% add dependence of the sensors 
    end

    preds=preds+0.05*coreP;    

    %% calculate sequenceness 
    for iShuf = 1:nShuf
        rp = uniquePerms(iShuf,:);  % use the 30 unique permutations (is.nShuf should be set to 29)
        T1 = TF(rp,rp); T2 = T1'; % backwards is transpose of forwards
        X=preds;
        
        nbins=maxLag+1;

       warning off
       dm=[toeplitz(X(:,1),[zeros(nbins,1)])];
       dm=dm(:,2:end);
       
       for kk=2:nstates
           temp=toeplitz(X(:,kk),[zeros(nbins,1)]);
           temp=temp(:,2:end);
           dm=[dm temp]; 
       end
       
       warning on
       
       Y=X;       
       betas = nan(nstates*maxLag, nstates);
%           
      %% GLM: state regression, with other lages       
       bins=maxLag;

       for ilag=1:bins
           temp_zinds = (1:bins:nstates*maxLag) + ilag - 1; 
           temp = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*Y;
           betas(temp_zinds,:)=temp(1:end-1,:);           
       end  

       betasnbins64=reshape(betas,[maxLag nstates^2]);
       bbb=pinv([T1(:) T2(:) squash(eye(nstates)) squash(ones(nstates))])*(betasnbins64'); %squash(ones(nstates))

       sf{iSj}(1,iShuf,2:end) = bbb(1,:); 
       sb{iSj}(1,iShuf,2:end) = bbb(2,:); 

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      %% Cross-Correlation
      for iLag=1:maxLag
          sf2{iSj}(1,iShuf,iLag+1) = sequenceness_Crosscorr(preds, T1, [], iLag);
          sb2{iSj}(1,iShuf,iLag+1) = sequenceness_Crosscorr(preds, T2, [], iLag);
      end        
                 
    end
end

sf = cell2mat(sf);
sb = cell2mat(sb);

sf2 = cell2mat(sf2);
sb2 = cell2mat(sb2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, 

%% GLM (fwd-bkw)
subplot(2,3,1)
npThresh = squeeze(max(abs(mean(sf(:,2:end,2:end)-sb(:,2:end,2:end),1)),[],3));
npThreshAll = max(npThresh);  
dtp = squeeze(sf(:,1,:)-sb(:,1,:));
shadedErrorBar(cTime, mean(dtp), std(dtp)/sqrt(nSubj)), hold on,
plot([cTime(1) cTime(end)], -npThreshAll*[1 1], 'k--'), plot([cTime(1) cTime(end)], npThreshAll*[1 1], 'k--')
title('GLM: fwd-bkw'), xlabel('lag (ms)'), ylabel('fwd minus bkw sequenceness')

%% GLM (fwd)
subplot(2,3,2)
npThresh = squeeze(max(abs(mean(sf(:,2:end,2:end),1)),[],3));
npThreshAll = max(npThresh);  
dtp = squeeze(sf(:,1,:));
shadedErrorBar(cTime, mean(dtp), std(dtp)/sqrt(nSubj)), hold on,
plot([cTime(1) cTime(end)], -npThreshAll*[1 1], 'k--'), plot([cTime(1) cTime(end)], npThreshAll*[1 1], 'k--')
title('GLM: fwd'), xlabel('lag (ms)'), ylabel('fwd sequenceness')

%% GLM (bkw)
subplot(2,3,3)
npThresh = squeeze(max(abs(mean(sb(:,2:end,2:end),1)),[],3));
npThreshAll = max(npThresh);  
dtp = squeeze(sb(:,1,:));
shadedErrorBar(cTime, mean(dtp), std(dtp)/sqrt(nSubj)), hold on,
plot([cTime(1) cTime(end)], -npThreshAll*[1 1], 'k--'), plot([cTime(1) cTime(end)], npThreshAll*[1 1], 'k--')
title('GLM: bkw'), xlabel('lag (ms)'), ylabel('bkw sequenceness')

%% Cross-Correlation (fwd-bkw)
sf=sf2;
sb=sb2;
subplot(2,3,4)
npThresh = squeeze(max(abs(mean(sf(:,2:end,2:end)-sb(:,2:end,2:end),1)),[],3));
npThreshAll = max(npThresh);  
dtp = squeeze(sf(:,1,:)-sb(:,1,:));
shadedErrorBar(cTime, mean(dtp), std(dtp)/sqrt(nSubj)), hold on,
plot([cTime(1) cTime(end)], -npThreshAll*[1 1], 'k--'), plot([cTime(1) cTime(end)], npThreshAll*[1 1], 'k--')
title('Correlation: fwd-bkw'), xlabel('lag (ms)'), ylabel('fwd minus bkw sequenceness')

%% Cross-Correlation (fwd)
subplot(2,3,5)
npThresh = squeeze(max(abs(mean(sf(:,2:end,2:end),1)),[],3));
npThreshAll = max(npThresh); 
dtp = squeeze(sf(:,1,:));
shadedErrorBar(cTime, mean(dtp), std(dtp)/sqrt(nSubj)), hold on,
plot([cTime(1) cTime(end)], -npThreshAll*[1 1], 'k--'), plot([cTime(1) cTime(end)], npThreshAll*[1 1], 'k--')
title('Correlation: fwd'), xlabel('lag (ms)'), ylabel('fwd sequenceness')

%% Cross-Correlation (bkw)
subplot(2,3,6)
npThresh = squeeze(max(abs(mean(sb(:,2:end,2:end),1)),[],3));
npThreshAll = max(npThresh); 
dtp = squeeze(sb(:,1,:));
shadedErrorBar(cTime, mean(dtp), std(dtp)/sqrt(nSubj)), hold on,
plot([cTime(1) cTime(end)], -npThreshAll*[1 1], 'k--'), plot([cTime(1) cTime(end)], npThreshAll*[1 1], 'k--')
title('Correlation: bkw'), xlabel('lag (ms)'), ylabel('bkw sequenceness')
