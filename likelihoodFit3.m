function [nlogL] = likelihoodFit3(sigma_m,sigma_p,sigma_s,expData,t,distTestEndpts,lookUpTab,maxDistAll,sigMmax,distFromTarget,endPtsFB,target1,reaches,indicated)

%PROSPECTIVE MODEL - uses only motor and setting noise

%All data converted to mm prior to being run through this model. All
%calculations assume inputs to be in mm.

%The inputs for this function are:
%sigma_m = test motor noise
%sigma_s = test setting noise
%trialNum = number of trials
%endpts = end points from data
%expData = circle sizes from data
%t = target location (centered data)
%matNan = maxtrix the size of the tablet in mm
%lookUPTab = look up table for MEG for this model
%maxDistAll = vector of all possible distances from target
%sigMmax = vector of possible sigma_m values
%distFromTarget = matrix the size of tablet in mm contianing distance from target at
%each mm

if nargin > 11
     %% CONTROL EXPERIMENT PARAMETER LIKELIHOOD
    
    RmTemp = 1/sigma_m^2;
    
    % log likelihood of sigma_m given target and endpoint:
    LLe = log2disotnormal(reaches,target1,sigma_m);
    
    RpTemp = 1/sigma_p^2;
    meanigivene = (RpTemp/(RpTemp+RmTemp))*reaches + ...
        (RmTemp/(RpTemp+RmTemp))*target1;
    SDigivene = (RpTemp/(RpTemp+RmTemp))*sigma_p;
    
    % log likelihood of sigma_p given sensed location and endpoint:
    LLigivene = log2disotnormal(indicated,meanigivene,SDigivene);
    
    % log likelihood of the sigma_p/sigma_m pair:
    LL = LLe + LLigivene;
end


%% INITIALIZE MEG CIRCLE SIZE
circSizeMEG = interp1(sigMmax,lookUpTab(1:length(sigMmax)),sigma_m);

%initalizing loop variables
probEndPt = nan(1,length(maxDistAll));
%% PROBABILITY OF ENDPOINTS FOR THESE PARAMETERS

%Probability of endpoint by distance from target
endPointProbDist = gauss2d(nan([268 476]), sigma_m, t);

for ii = 1:length(maxDistAll)                       %one through maximum distance
    probEndPt(ii) = sum(endPointProbDist(round(distFromTarget) == ii));    %probabilites in the signal distribution for a given distance.
end

circDist = normpdf(maxDistAll,circSizeMEG,sigma_s); %find the distribution around the MEG circle using sigma_s

%% PROBABILITIES OF CIRCLES

probTestEndpts = interp1(maxDistAll,probEndPt,distTestEndpts);       %probabilites of endpoints at a given distance
probFBEndpts = interp1(maxDistAll,probEndPt,endPtsFB);

circleProb = interp1(maxDistAll,circDist,expData);                %compare the test circle size to the distribution around the MEG circle

%% COMBINE PROBABILITES
nLL = -(sum(log(probTestEndpts)) + sum(log(probFBEndpts)) + sum(log(circleProb))); %add the log likelihood of the end points to the log likelihood of the circles

if nargin >11
    nlogL = nLL - LL;
else
    nlogL = nLL;
end
%% Helper Functions
    function matReturn = gauss2d(mat, sigma, center)
        %Takes 'mat', a matrix of the size of the tablet space in mm 
        %Values in this 'mat' matrix do not matter. 'sigma' is the SD of distribution
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
