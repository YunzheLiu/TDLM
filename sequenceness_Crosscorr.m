function [sf] = sequenceness_Crosscorr(rd, T, T2, lag)
% rd is samples by states 
% T is the transition matrix of interest
% T2 is the 2-step transition matrix, or can be []
% lag is how many samples the data should be shifted by 

if isempty(T2)
    nSamples = size(rd,1); nStates = size(rd,2);
    orig = rd(1:(end-2*lag),:)*T;
    proj = rd((1+lag):(end-lag),:);

   %% Scale variance 
   corrtemp=nan(size(proj,2),1) ;
   for iseq=1:size(proj,2)
       if nansum(orig(:,iseq))~=0 && nansum(proj(:,iseq))~=0
           corrtemp(iseq)=corr(orig(:,iseq),proj(:,iseq));
       end
   end
   sf = nanmean(corrtemp);
    
else
    nSamples = size(rd,1); nStates = size(rd,2);
    orig = rd(1:(end-2*lag),:)*T2;
    proj = rd((1+lag):(end-lag),:)*T;
    proj2 = rd((1+2*lag):end,:);
    sf = nanmean(nansum((orig - nanmean(orig,2)) .* (proj - nanmean(proj,2)) .* (proj2 - nanmean(proj2,2)),1));
end

end




