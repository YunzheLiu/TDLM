function Pr = placeBayes(X, rateMaps, binsize)

% INPUT
% X - candiate event to be decoded, ntime * ncells, each entry is the spike counts
% rateMaps - rate maps, ncells * npos
% binSize, in the unit of s, this is tau in the equation

% OUTPUT
% Pr, decoded posterior

Pr=[];
expterm = exp(-binsize*sum(rateMaps,1)); % sum over neurons
prior   = 1/size(rateMaps,2); % uniform over space

for ipos=1:size(rateMaps,2)

    for iT=1:size(X,1)
        temp=[];
        for icell=1:size(rateMaps,1)
            temp(icell)=rateMaps(icell,ipos)^X(iT,icell);
        end
        
        Pr(iT,ipos) = prior*prod(temp)*expterm(ipos);
    end
end

end