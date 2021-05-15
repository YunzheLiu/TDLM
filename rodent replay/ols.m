function [cope,varcope,tstat]=ols(data,des,tc)
% [COPE,VARCOPE]=OLS(DATA,DES,TC)
% DATA IS T x V
% DES IS T x EV (design matrix)
% TC IS NCONTRASTS x EV  (contrast matrix)
%
% TB 2004
if (nargin<3)
    tc=eye(size(des,2));
end

if(size(data,1)~=size(des,1))
  error('OLS::DATA and DES have different number of time points');
elseif(size(des,2)~=size(tc,2))
  error('OLS:: DES and EV have different number of evs')
end


%remove NaNs and warn
test=sum(data,2)+sum(des,2);
if(isnan(sum(test))); warning('Removing NaN rows');
des(isnan(test),:)=[]; data(isnan(test),:)=[];
end


pdes=pinv(des);
prevar=diag(tc*pdes*pdes'*tc');


% R=eye(size(des,1))-des*pdes;
% tR=trace(R);
tR=size(des,1)-size(des,2);

pe=pdes*data;
cope=tc*pe; 
if(nargout>1)
  res=data-des*pe;
  sigsq=sum(res.*res/tR);       
 
  varcope=prevar*sigsq;
  
  if(nargout>2)
    tstat=cope./sqrt(varcope);
  end
end


