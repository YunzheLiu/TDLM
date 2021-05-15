function [Msf, Msb, Mpf, Mpb, Tsf, Tsb, Tpf, Tpb] = TDLM_Pweight_wTime(X,T,nstates,nlags,nshuf)

% in this version, we also model time as well interaction with it in the
% first level GLM

timeR = normalize(1:size(X,1));

[~, uniquePerms] = uperms([1:nstates],nshuf);
nbins=nlags+1;

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

Msf=nan(1,nshuf,nlags);
Msb=nan(1,nshuf,nlags);
Mpf=nan(1,nshuf,nlags);
Mpb=nan(1,nshuf,nlags);

Tsf=nan(1,nshuf,nlags);
Tsb=nan(1,nshuf,nlags);
Tpf=nan(1,nshuf,nlags);
Tpb=nan(1,nshuf,nlags);

if size(uniquePerms,1)<nshuf
    nshuf=size(uniquePerms,1);
end
    
%% GLM: state regression, with other lages 
betas1 = nan(nstates*nlags, nstates);
vars1  = nan(nstates*nlags, nstates);

betas2 = nan(nstates*nlags, nstates);
vars2  = nan(nstates*nlags, nstates);

bins  = nlags;

for ilag=1:bins
    temp_zinds = (1:bins:nstates*nlags) + ilag - 1; 
    cdesgin=dm(:,temp_zinds);   
    
    %% modelling time
    XR=[];
    for istate=1:size(cdesgin,2)
        XR(:,istate) = cdesgin(:,istate).*timeR';
    end
    
    %% GLM
    for istate=1:size(Y,2)

        [cope,varcope,~]=ols(Y(:,istate),[ones(length(cdesgin),1),cdesgin, XR]);

        betas1(temp_zinds,istate) = cope(2:1+size(cdesgin,2));
        vars1(temp_zinds,istate)  = varcope(2:1+size(cdesgin,2));
        
        betas2(temp_zinds,istate) = cope(2+size(cdesgin,2):end);
        vars2(temp_zinds,istate)  = varcope(2+size(cdesgin,2):end);  
    end
        
end

betaM      = reshape(betas1,[nlags nstates^2]);
precisionM = 1./reshape(vars1,[nlags nstates^2]);
weightedM  = betaM.*precisionM;

betaT      = reshape(betas2,[nlags nstates^2]);
precisionT = 1./reshape(vars2,[nlags nstates^2]);
weightedT  = betaT.*precisionT;

Tall= squash(ones(nstates)) - squash(eye(nstates));

for iShuf = 1:nshuf
    rp = uniquePerms(iShuf,:);
    T1 = T(rp,rp); 
    T2 = T1'; % backwards is transpose of forwards 
    
    for ilag=1:nlags
       
        %% sequence effect 

        Sf_M  = nansum(weightedM(ilag,:)'.*T1(:))/nansum(precisionM(ilag,:)'.*T1(:));
        Sb_M  = nansum(weightedM(ilag,:)'.*T2(:))/nansum(precisionM(ilag,:)'.*T2(:)); 
        Smean = nansum(weightedM(ilag,:)'.*Tall)/nansum(precisionM(ilag,:)'.*Tall);

        Msf(1,iShuf,ilag) = Sf_M-Smean;
        Msb(1,iShuf,ilag) = Sb_M-Smean; 
        
        Mpf(1,iShuf,ilag) = nansum(precisionM(ilag,:)'.*T1(:)); 
        Mpb(1,iShuf,ilag) = nansum(precisionM(ilag,:)'.*T2(:));          
        
        %% time modulation effect

        Sf_T  = nansum(weightedT(ilag,:)'.*T1(:))/nansum(precisionT(ilag,:)'.*T1(:));
        Sb_T  = nansum(weightedT(ilag,:)'.*T2(:))/nansum(precisionT(ilag,:)'.*T2(:)); 
        Tmean = nansum(weightedT(ilag,:)'.*Tall)/nansum(precisionT(ilag,:)'.*Tall);

        Tsf(1,iShuf,ilag) = Sf_T-Tmean;
        Tsb(1,iShuf,ilag) = Sb_T-Tmean; 
        
        Tpf(1,iShuf,ilag) = nansum(precisionT(ilag,:)'.*T1(:)); 
        Tpb(1,iShuf,ilag) = nansum(precisionT(ilag,:)'.*T2(:));   
        
    end
end
