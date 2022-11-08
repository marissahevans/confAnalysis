function [nlogL] = likelihoodFitRetro(sigma_m,sigma_p,sigma_s,trialNum,endpts,expData,t,distTestEndpts,lookUpTab,maxDistAll,sigPmax,distFromTarget,endPtsFB,target1,reaches,indicated)

%RETROSPECTIVE MODEL - uses motor, prosprioceptive, and setting noise, but
%motor noise is not used after the endpoints are established

%All data converted to mm prior to being run through this model. All
%calculations assume inputs to be in mm.

%The inputs for this function are:
%sigma_m = test motor noise
%sigma_p = test proprioceptive noise
%sigma_s = test setting noise
%trialNum = number of trials
%endpts = end points from data
%simData = circle sizes from data
%t = target location (for centered data)
%matNan = maxtrix the size of the tablet in mm
%lookUPTab = look up table for MEG for this model
%maxDistAll = vector of all possible distances from target
%sigPmax = vector of possible sigma_p values
%distFromTarget = matrix the size of tablet in mm contianing distance from target at
%each mm

if nargin > 13
    %% CONTROL EXPERIMENT PARAMETER LIKELIHOOD
    
    % log likelihood of sigma_m given target and endpoint:
    LLe = log2disotnormal(reaches,target1,sigma_m);
        
    % log likelihood of sigma_p given sensed location and endpoint:
    LLigivene = log2disotnormal(indicated,reaches,sigma_p);
    
    % log likelihood of the sigma_p/sigma_m pair:
    LL = LLe + LLigivene;
end

circSizeMEG = nan(1,length(lookUpTab));
%% INITIALIZE MEG CIRCLE SIZE
for ii = 1:length(lookUpTab)
    circSizeMEG(ii) = interp1(sigPmax,lookUpTab(ii,1:length(sigPmax)),sigma_p);
end

%initalizing loop variables
probEndPt = nan(1,length(maxDistAll));
sumCirc = nan(1,trialNum);
circleProb = nan([268 476]);
%% PROBABILITY OF ENDPOINTS FOR THESE PARAMETERS

%log likelihood of endpoint by distance from target
endPointProbDist = gauss2d(nan([268 476]), sigma_m, t);

for ii = 1:length(maxDistAll)                       %one through maximum distance
    probEndPt(ii) = sum(endPointProbDist(round(distFromTarget) == ii));    %probabilites in the signal distribution for a given distance.
end

probTestEndpts = interp1(maxDistAll,probEndPt,distTestEndpts);       %probabilites of endpoints at a given distance
probFBEndpts = interp1(maxDistAll,probEndPt,endPtsFB);

%% Step 2 - FIND DISTANCES FROM EACH P TO TARGET
for tt = 1:trialNum %begin looping over all trials
    
    %End point on this trial from data
    e(1) = endpts(tt,1);
    e(2) = endpts(tt,2);
    
    circ = expData(tt);
    
    %log likelihood distribution of sensed locations
    signalDist = gauss2d(nan([268 476]), sigma_p, e);
    
    
    %% Step 3 - PROBABILITIES OF CIRCLES
    %Calculate probability of a circle size being selected at a given
    %signal distance
    
    for ii = 1:size(distFromTarget,1)               %loop through all distances x
        for jj = 1:size(distFromTarget,2)           %loop through all distances y
            circleOpt = interp1(maxDistAll,circSizeMEG,distFromTarget(ii,jj));%find MEG circle for this specific distance
            circleProb(ii,jj) = normpdf(circ,circleOpt,sigma_s);     %compare the circle from the data to the probability of the MEG circle
        end
    end
    
    %% STEP 4 - COMBINE PROBABILITES
    combinedDist = circleProb .* signalDist;              %combine the circle probability with the signal probability
    
    sumCirc(tt) = sum(combinedDist(:));                   %sum matrix for each trial
    
end %trial loop
nLL = -(sum(log(probTestEndpts)) + sum(log(probFBEndpts)) + sum(log(sumCirc))); %add the log likelihood of the end points to the log likelihood of the circles

if nargin >13
    nlogL = nLL - LL;
else
    nlogL = nLL;
end
%% Helper Functions
    function matReturn = gauss2d(mat, sigma, center)
        %Takes 'mat', a matrix of the size you want the final output space to be
        %in. Values in this matrix do not matter. 'sigma' is the SD of distribution
        %and 'center' is the peak of the returned gaussian (internal function
        %gaussC is necessary to run this function).
        
        gsize = size(mat);                      %find size of sample space
        [R,C] = ndgrid(1:gsize(1), 1:gsize(2)); %create a grid of the sample space size
        matReturn = gaussC(C,R, sigma, center); %function output with imbedded function gaussC
    end

    function val = gaussC(x, y, sigma, center)
        %Returns x by y matrix with a gaussian distribution with peak at center,
        %and SD of sigma.
        
        xc = center(1);                         %x coordinate of centr
        yc = center(2);                         %y coordinate of center
        exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma^2); %gaussian exponent
        amplitude = 1 / (sigma^2 * 2*pi);       %constant
        val       = amplitude * exp(-exponent); %gaussian 2d
    end

% direct computation of log(p(data|mu,sigma))
% data is Nx2
% mu is 1x2 or Nx2
% sigma is a scalar SD of the isotropic 2-d Gaussian

    function ll = log2disotnormal(data,mu,sigma)
        
        centered = data - mu;
        N = size(data,1);
        ll = -N*log(2*pi) - 2*N*log(sigma) - (sum(centered(:).^2)/(2*sigma^2));
    end
end

