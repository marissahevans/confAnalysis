%% fmincon for retrospective model #3

load simOutput
load model3LookUpTab

%Tablet specs
radiusWac = linspace(0,69,150);

xCen = 238;                                 %center of tablet x coordinate
yCen = 134;                                 %cetner of tablet y coordinate
t = [xCen yCen];                            %target location for all trials
tabSize = [268 476];                        %outside bounds of tablet space
matNan = nan(tabSize);                      %matrix of nans the size of the tablet
sigMmax = 1:100;                           %vector of test sigma_m values

% Find all distances to target on tablet 
for nn = 1:size(matNan,1)
    for mm = 1:size(matNan,2)
        distFromTarget(nn,mm) = sqrt((mm - t(1))^2 + (nn-t(2))^2);
    end
end

modelFitCircleNoise1(modelFitCircleNoise1 <1) = 1;
modelFitCircleNoise2(modelFitCircleNoise2 <1) = 1;
modelFitCircleNoise3(modelFitCircleNoise3 <1) = 1;

maxDistAll = 0:max(distFromTarget(:)); %range of distance to max distance on the tablet. 
distFromTarget(distFromTarget > max(maxDistAll)) = max(maxDistAll);

%fmincon parameters
lb = [1, 1];                          %lower bound
ub = [40,10];                         %upper bound
init = rand*(ub-lb)+lb;                  %initiation
options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off'); %fmincon options

numTrials = [100,200,300,400,500,600,700,800];
numSims = 20;

for jj = 1:length(numTrials)
    
    for ii = 1:numSims
        distTestEndpts1 = sqrt((ePts1X(ii,1:numTrials(jj)) - t(1)).^2 + (ePts1Y(ii,1:numTrials(jj))-t(2)).^2); %distances from endpoints to target
        distTestEndpts2 = sqrt((ePts2X(ii,1:numTrials(jj)) - t(1)).^2 + (ePts2Y(ii,1:numTrials(jj))-t(2)).^2); %distances from endpoints to target
        distTestEndpts3 = sqrt((ePts3X(ii,1:numTrials(jj)) - t(1)).^2 + (ePts3Y(ii,1:numTrials(jj))-t(2)).^2); %distances from endpoints to target
        
        m3nLogLA = @(p) likelihoodFit3(p(1),[],p(2),modelFitCircleNoise1(ii,1:numTrials(jj)),t,distTestEndpts1,fit3LookUpMat,maxDistAll,sigMmax,distFromTarget);
        m3nLogLB = @(p) likelihoodFit3(p(1),[],p(2),modelFitCircleNoise2(ii,1:numTrials(jj)),t,distTestEndpts2,fit3LookUpMat,maxDistAll,sigMmax,distFromTarget);
        m3nLogLC = @(p) likelihoodFit3(p(1),[],p(2),modelFitCircleNoise3(ii,1:numTrials(jj)),t,distTestEndpts3,fit3LookUpMat,maxDistAll,sigMmax,distFromTarget);
        
        [estP3A(ii,:), minNLL3A(ii)] = fmincon(m3nLogLA,init,[],[],[],[],lb,ub,[],options); %find minimums - model 3 w/ data 1
        [estP3B(ii,:), minNLL3B(ii)] = fmincon(m3nLogLB,init,[],[],[],[],lb,ub,[],options); %find minimums - model 3 w/ data 2
        [estP3C(ii,:), minNLL3C(ii)] = fmincon(m3nLogLC,init,[],[],[],[],lb,ub,[],options); %find minimums - model 3 w/ data 3
        
    end
    sigMest1(:,jj) = estP3A(:,1);
    sigMest2(:,jj) = estP3B(:,1);
    sigMest3(:,jj) = estP3C(:,1);
    
    sigSest1(:,jj) = estP3A(:,2);
    sigSest2(:,jj) = estP3B(:,2);
    sigSest3(:,jj) = estP3C(:,2);
    jj
    
end


save model3output.mat sigMest1 sigMest2 sigMest3 sigSest1 sigSest2 sigSest3 minNLL3A minNLL3B minNLL3C
