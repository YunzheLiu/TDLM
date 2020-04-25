clear;
clc;
close all;

%%  transition matrix
TF = [0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0];
TR = TF';
nSensors = 273;
nTrainPerStim = 18; % how many training examples (direct presentations) for each stimulus)
nNullExamples = nTrainPerStim*8; % how many null examples to use
nSamples = 8000; % 60 seconds of unlabelled data to predict
nSequences = 6000; % how many real sequences to put in the data
maxLag = 60; % evaluate cross-correlation up to 600ms
cTime = 0:10:maxLag*10; % the milliseconds of each cross-correlation time lag
nSubj = 24; % number of subjects to simulate
gamA = 10;%10;%1e4; 
gamB = 0.5;%0.6;%;6e-4; % parameters for the gamma distribution of intervals between states in a sequence
nstates=8;
[~, pInds] = uperms([1:8],29);
uniquePerms=pInds;
nShuf = size(uniquePerms,1);
samplerate=100;
nlength=3; % change here to test for multi-length

sf1 = cell(nSubj,1);  sb1 = cell(nSubj,1);
sf2 = cell(nSubj,1);  sb2 = cell(nSubj,1);  sc2 = cell(nSubj,1);

parfor iSj = 1:nSubj

    sf1{iSj} = nan(1, nShuf, maxLag+1);
    sb1{iSj} = nan(1, nShuf, maxLag+1);   
    
    sf2{iSj} = nan(1, nShuf, maxLag+1);
    sb2{iSj} = nan(1, nShuf, maxLag+1);  
    sc2{iSj} = nan(1, nShuf, maxLag+1);
    
    disp(['iSj=' num2str(iSj)])
    
   %% generate dependence of the sensors
    A = randn(nSensors);
    [U,~] = eig((A+A')/2); 
    covMat = U*diag(abs(randn(nSensors,1)))*U';
    
    %% generate the true patterns
    commonPattern = randn(1,nSensors);    
    patterns =  repmat(commonPattern, [8 1]) + randn(8, nSensors); % this induce correlations between classifers, for seperate representation: randn(8, nSensors);
   
    %% make training data
    trainingData = 4*randn(nNullExamples+8*nTrainPerStim, nSensors) + [zeros(nNullExamples,nSensors); ...
        repmat(patterns, [nTrainPerStim 1])];  
    trainingLabels = [zeros(nNullExamples,1); repmat((1:8)', [nTrainPerStim 1])];

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
        state = false(8,1); 
%         state(randsample([4,8],1)) = true;  % start the sequence in a random state
        state(randsample([1,5],1)) = true;  % start the sequence in a random state

        for iMv=1:nlength
            if sum(state)~=0
                X(seqTime,:) = X(seqTime,:) + patterns(state,:);%randn(1,nSensors); %cos(pi./randsample([7:13],nSensors,'true'));%patterns(state,:);%+ cos(pi./randsample([7:13],nSensors,'true'));
                state = (state'*TF)'; state2 = false(8,1); state2(find(rand < cumsum(state), 1, 'first')) = true; state = state2; % advance states
                seqTime = seqTime + round(gamrnd(gamA,gamB));
            end
        end
    end

    %% make predictions with trained models
     preds = 1./(1+exp(-(X*betas + repmat(intercepts, [nSamples 1]))));
    
    %% Prepare T matrix    
    Tf_Y =[3,4,7,8];    
    Tf_X2=[2,3,6,7];
    Tf_X1=[1,2,5,6];
    
    TF2=zeros(length(Tf_Y),nstates);
    TFauto=zeros(length(Tf_Y),nstates);
    for ist=1:length(Tf_X1)
        TF2(ist,Tf_Y(ist))=1;
        TFauto(ist,unique([Tf_X1(ist),Tf_X2(ist)]))=1;
    end
        
    Tr_Y =[1,2,5,6];
    Tr_X2=[2,3,6,7];
    Tr_X1=[3,4,7,8];
    TR2=zeros(length(Tr_Y),nstates);
    TRauto=zeros(length(Tr_Y),nstates);
    for ist=1:length(Tr_X1)
        TR2(ist,Tr_Y(ist))=1;
        TRauto(ist,unique([Tr_X1(ist),Tr_X2(ist)]))=1;
    end    
           
    %% Core sequence detection
    X=preds;
    Y=X;

    X2bin=nan(maxLag, nSamples, nstates, nstates);  

    for ilag=1:maxLag
        pad = zeros([ilag nstates]);         
        X1 = [pad; pad; X(1:end-2*ilag,:)];
        X2 = [pad; X(1:end-ilag,:)];

        for i=1:nstates
            for j=1:nstates
                X2bin(ilag,:,i,j)=X1(:,i).*X2(:,j);
            end
        end
    end

    betaF=nan(maxLag,length(Tf_Y),nstates);
    betaR=nan(maxLag,length(Tr_Y),nstates);

    for ilag=1:maxLag        
        Xmatrix=squeeze(X2bin(ilag,:,:,:));

        for istate=1:length(Tf_Y)
            Xfwd=squeeze(Xmatrix(:,:,Tf_X2(istate))); 
            Xbkw=squeeze(Xmatrix(:,:,Tr_X2(istate)));
% 
%             Xfwd(:,Tf_Y(istate))=0;
%             Xbkw(:,Tr_Y(istate))=0;

            tempF = pinv([Xfwd ones(length(Xfwd),1)])*Y;                
            betaF(ilag,istate,:)=tempF(Tf_X1(istate),:);  

            tempR = pinv([Xbkw ones(length(Xbkw),1)])*Y;                  
            betaR(ilag,istate,:)=tempR(Tr_X1(istate),:);             
        end
    end
    
    betaF=reshape(betaF,[maxLag,length(Tf_Y)*nstates]);
    betaR=reshape(betaR,[maxLag,length(Tr_Y)*nstates]);
        
    for iShuf = 1:nShuf
        rp = uniquePerms(iShuf,:); 
       
        %% 2nd level         
        constF=ones(length(Tf_Y),nstates);
        constR=ones(length(Tr_Y),nstates);      
        
        TFs=TF2(:,rp);
        TRs=TR2(:,rp); 
        
        cc=pinv([squash(TFs) squash(TFauto) squash(constF)])*betaF';     
        sf1{iSj}(1,iShuf,2:end)=cc(1,:); 
                       
        cc=pinv([squash(TRs) squash(TRauto) squash(constR)])*betaR';    
        sb1{iSj}(1,iShuf,2:end)=cc(1,:);
        
    end    
end

sf1 = cell2mat(sf1);
sb1 = cell2mat(sb1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure, 

%% calculate threshold and plot (2 fwd)
subplot(1,2,1)
npThresh = squeeze(max(abs(mean(sf1(:,2:end,2:end),1)),[],3));
npThreshAll = max(npThresh);  % can also use max instead of 95%, very similar
dtp = squeeze(sf1(:,1,:));
shadedErrorBar(cTime(1:end), mean(dtp), std(dtp)/sqrt(nSubj)), hold on,
plot([cTime(1) cTime(end)], -npThreshAll*[1 1], 'k--'), plot([cTime(1) cTime(end)], npThreshAll*[1 1], 'k--')
title('fwd'), xlabel('lag (ms)'), ylabel('L2 - fwd sequenceness')

subplot(1,2,2)
npThresh = squeeze(max(abs(mean(sb1(:,2:end,2:end),1)),[],3));
npThreshAll = max(npThresh);  % can also use max instead of 95%, very similar
dtp = squeeze(sb1(:,1,:));
shadedErrorBar(cTime(1:end), mean(dtp), std(dtp)/sqrt(nSubj)), hold on,
plot([cTime(1) cTime(end)], -npThreshAll*[1 1], 'k--'), plot([cTime(1) cTime(end)], npThreshAll*[1 1], 'k--')
title('bkw'), xlabel('lag (ms)'), ylabel('L2 - bkw sequenceness')
