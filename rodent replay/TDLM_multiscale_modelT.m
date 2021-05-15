function [Seq_integrate, Seq_separate, Time_integrate, Time_separate]=TDLM_multiscale_modelT(X, rateMaps, binSize, nshuf, nlags, TOI)

% INPUT
% X          - candinate replay evernts, ntime (concatenated) * ncells
% rateMaps   - raw rate maps, ncells * npos (this inlcudes both inbound and outbound)  
% binSize    - temporal bin size, in s, so 0.001 -> 1ms
% nShuf      - number of nShufs
% nlags      - number of time lags to test
% TOI        - time lag of interest

% OUTPUT
% interSeq   - integrate sequence speed
% separaSeq  - save sequences strength in different speeds

ratemap=[];

%% scale and corrsponding time lags to average (hard coded for now)
scales=[24,12,8,6]; % numnber of states, higher number indciates higher spatial resolution
ntimes={TOI;TOI*2;TOI*3;TOI*4}; % corresponding time lags of interest, this is for 2 ms, transiting 2 cm, i.e., 1000cm/s

%% make the rate map in different scales - for marginalize over in-, out-bound
PosRaw=size(rateMaps,2)/2;

RawIn=rateMaps(:,1:PosRaw);
RawOut=rateMaps(:,PosRaw+1:end);

for iscale=1:length(scales)

    ratemap1=[];
    ratemap2=[];
    step=floor(PosRaw/scales(iscale));    
    
    for istate=1:scales(iscale)
        if istate<scales(iscale)
            stateRange=1+(istate-1)*step:istate*step;
        else
            stateRange=1+(istate-1)*step:PosRaw;
        end
        ratemap1(:,istate)=nanmean(RawIn(:,stateRange),2);
        ratemap2(:,istate)=nanmean(RawOut(:,stateRange),2);
    end    
    
    ratemap{iscale,1}=[ratemap1]; % inbound
    ratemap{iscale,2}=[ratemap2]; % outbound
    
end
     
%% Core of multiscale TDLM

Seq_separate  = [];
Seq_integrate = [];

Time_separate  = [];
Time_integrate = [];

Msfall=[];
Msball=[];
Mpfall=[];
Mpball=[];

Tsfall=[];
Tsball=[];
Tpfall=[];
Tpball=[];

for irun=1:2
    
    for iscale=1:size(ratemap,1)

        %% separately for different spatial scales
        Pr  = placeBayes(X, ratemap{iscale, irun}, binSize);
        Pr2 = Pr./(sum(Pr,2));  
        
        [rows, ~] = find(isnan(Pr2));
        Pr2(unique(rows),:)=[];
        
        %% sequence analysis core
        numPlace=size(Pr2,2);
        T=diag(ones(numPlace-1,1),1);
        
        if irun==2
            T=T'; % because I code up linearized state as inbound running
        end
        
        [Msf, Msb, Mpf, Mpb, Tsf, Tsb, Tpf, Tpb] = TDLM_Pweight_wTime(Pr2,T,numPlace,nlags,nshuf+1); % model time and its interaction explicitly in the first level

        %% Save sequence results in each scale separately
        Seq_separate (1,irun,iscale,:,:)=squeeze(Msf(1,:,:));
        Seq_separate (2,irun,iscale,:,:)=squeeze(Msb(1,:,:));
        
        Time_separate(1,irun,iscale,:,:)=squeeze(Tsf(1,:,:));
        Time_separate(2,irun,iscale,:,:)=squeeze(Tsb(1,:,:));
                
        Msfall(irun,iscale,:,1)=squeeze(Msf(1,:,ntimes{iscale}));
        Msball(irun,iscale,:,1)=squeeze(Msb(1,:,ntimes{iscale}));
        Mpfall(irun,iscale,:,1)=squeeze(Mpf(1,:,ntimes{iscale}));
        Mpball(irun,iscale,:,1)=squeeze(Mpb(1,:,ntimes{iscale})); 
        
        Tsfall(irun,iscale,:,1)=squeeze(Tsf(1,:,ntimes{iscale}));
        Tsball(irun,iscale,:,1)=squeeze(Tsb(1,:,ntimes{iscale}));
        Tpfall(irun,iscale,:,1)=squeeze(Tpf(1,:,ntimes{iscale}));
        Tpball(irun,iscale,:,1)=squeeze(Tpb(1,:,ntimes{iscale}));   
    end
    
    %% weighted average across scale     
    for ishuf=1:nshuf+1
        for ilag=1:size(Msfall,4)
            Seq_integrate(1, irun, ishuf,ilag)=nansum(Msfall(irun,:,ishuf,ilag).*Mpfall(irun,:,ishuf,ilag))/nansum(Mpfall(irun,:,ishuf,ilag));
            Seq_integrate(2, irun, ishuf,ilag)=nansum(Msball(irun,:,ishuf,ilag).*Mpball(irun,:,ishuf,ilag))/nansum(Mpball(irun,:,ishuf,ilag));

            Time_integrate(1, irun, ishuf,ilag)=nansum(Tsfall(irun,:,ishuf,ilag).*Tpfall(irun,:,ishuf,ilag))/nansum(Tpfall(irun,:,ishuf,ilag));
            Time_integrate(2, irun, ishuf,ilag)=nansum(Tsball(irun,:,ishuf,ilag).*Tpball(irun,:,ishuf,ilag))/nansum(Tpball(irun,:,ishuf,ilag));                
        end
    end
    disp(['finish irun=',num2str(irun)])
  
end

end