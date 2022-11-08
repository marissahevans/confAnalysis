function [sigPmarg,sigMmarg] = ctrlFit(target,reaches,indicated,m_vec,p_vec)

% Simultaneously estimate sigma_m and sigma_p by ML
LL=zeros(length(m_vec),length(p_vec));
LLigivene=zeros(length(m_vec),length(p_vec));
LLe=zeros(length(m_vec),1);

for ii = 1:length(m_vec)        %loop over all sigma_m options
    RmTemp = 1/m_vec(ii)^2;
    
    % log likelihood of sigma_m given target and endpoint:
    LLe(ii) = log2disotnormal(reaches,target,m_vec(ii));
    
    for jj = 1:length(p_vec)    %loop over all sigma_pß options
        RpTemp = 1/p_vec(jj)^2;
        meanigivene = (RpTemp/(RpTemp+RmTemp))*reaches + ...
            (RmTemp/(RpTemp+RmTemp))*target;
        SDigivene = (RpTemp/(RpTemp+RmTemp))*p_vec(jj);
        % log likelihood of sigma_p given sensed location and endpoint:
        LLigivene(ii,jj) = log2disotnormal(indicated,meanigivene,SDigivene);
        % log likelihood of the sigma_p/sigma_m pair:
        LL(ii,jj) = LLe(ii) + LLigivene(ii,jj);
    end
end

NormPost = exp(LL - max(LL(:)));

mmarg = sum(NormPost,2);
mmarg = mmarg/sum(mmarg);
sigMmarg = m_vec(find(mmarg == max(mmarg)));

pmarg = sum(NormPost,1);
pmarg = pmarg/sum(pmarg);
sigPmarg = p_vec(find(pmarg == max(pmarg)));


    function ll = log2disotnormal(data,mu,sigma)
        
        centered = data - mu;
        N = size(data,1);
        ll = -N*log(2*pi) - 2*N*log(sigma) - (sum(centered(:).^2)/(2*sigma^2));
    end
end
