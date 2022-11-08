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
sigPmax = 1:100;                           %vector of test sigma_p values


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
lb = 1;                          %lower bound
ub = 50;                         %upper bound
init = rand*(ub-lb)+lb;                  %initiation
options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off'); %fmincon options

numTrials = 300;
numSims = 20;



for ii = 1:numSims
    [sigP3A(ii),sigM3A(ii)] = ctrlFit(target1,[reach1X(ii,:); reach1Y(ii,:)]',[indic1X(ii,:); indic1Y(ii,:)]',sigMmax,sigPmax);
    [sigP3B(ii),sigM3B(ii)] = ctrlFit(target1,[reach2X(ii,:); reach2Y(ii,:)]',[indic2X(ii,:); indic2Y(ii,:)]',sigMmax,sigPmax);
    [sigP3C(ii),sigM3C(ii)] = ctrlFit(target1,[reach3X(ii,:); reach3Y(ii,:)]',[indic3X(ii,:); indic3Y(ii,:)]',sigMmax,sigPmax);
    
    distTestEndpts1 = sqrt((ePts1X(ii,1:numTrials) - t(1)).^2 + (ePts1Y(ii,1:numTrials)-t(2)).^2); %distances from endpoints to target
    distTestEndpts2 = sqrt((ePts2X(ii,1:numTrials) - t(1)).^2 + (ePts2Y(ii,1:numTrials)-t(2)).^2); %distances from endpoints to target
    distTestEndpts3 = sqrt((ePts3X(ii,1:numTrials) - t(1)).^2 + (ePts3Y(ii,1:numTrials)-t(2)).^2); %distances from endpoints to target
    
    m3nLogLA = @(p) likelihoodFit3(sigM3A(ii),sigP3A(ii),p,modelFitCircleNoise1(ii,1:numTrials),t,distTestEndpts1,fit3LookUpMat,maxDistAll,sigMmax,distFromTarget);
    m3nLogLB = @(p) likelihoodFit3(sigM3B(ii),sigP3B(ii),p,modelFitCircleNoise2(ii,1:numTrials),t,distTestEndpts2,fit3LookUpMat,maxDistAll,sigMmax,distFromTarget);
    m3nLogLC = @(p) likelihoodFit3(sigM3C(ii),sigP3C(ii),p,modelFitCircleNoise3(ii,1:numTrials),t,distTestEndpts3,fit3LookUpMat,maxDistAll,sigMmax,distFromTarget);
    
    [sigS3A(ii), minNLL3A(ii)] = fmincon(m3nLogLA,init,[],[],[],[],lb,ub,[],options); %find minimums - model 3 w/ data 1
    [sigS3B(ii), minNLL3B(ii)] = fmincon(m3nLogLB,init,[],[],[],[],lb,ub,[],options); %find minimums - model 3 w/ data 2
    [sigS3C(ii), minNLL3C(ii)] = fmincon(m3nLogLC,init,[],[],[],[],lb,ub,[],options); %find minimums - model 3 w/ data 3
    
end

estP3A = [sigM3A', sigP3A', sigS3A'];
estP3B = [sigM3B', sigP3B', sigS3B'];
estP3C = [sigM3C', sigP3C', sigS3C'];

save model3output_fixP.mat estP3A estP3B estP3C minNLL3A minNLL3B minNLL3C
